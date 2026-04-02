import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import * as ngl from 'NGL';

import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {IPdbHelper, PdbResDataFrameType} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {Molecule3DUnits} from '@datagrok-libraries/bio/src/molecule-3d/molecule-3d-units-handler';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {AtomBase, AtomCoordsBase, LineBase} from '@datagrok-libraries/bio/src/pdb/format/types-base';
import {PdbAtomCoords, PdbAtomTer} from '@datagrok-libraries/bio/src/pdb/format/types-pdb';

import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {DefaultPluginSpec, PluginSpec} from 'molstar/lib/mol-plugin/spec';
import {PluginConfig} from 'molstar/lib/mol-plugin/config';
import {StateObjectSelector} from 'molstar/lib/mol-state';
import {Model} from 'molstar/lib/mol-model/structure';
import {PluginStateObject} from 'molstar/lib/mol-plugin-state/objects';
import {Sequence} from 'molstar/lib/mol-model/sequence';
import {Pdbqt} from './pdbqt-parser';
import {IMPORT} from '../consts-import';


// -- PDB header parsing types and functions --
// Parses metadata from PDB file headers (TITLE, SOURCE, JRNL, SITE, etc.)
// Complementary to Molstar which handles coordinates/sequences but not header metadata.
// All fields are optional to handle files from RCSB, MOE, Schrödinger, and other sources.

export interface PdbEntity {
  molId: string;
  molecule?: string;
  chains?: string[];
  engineered?: boolean;
  mutation?: boolean;
}

export interface PdbMutation {
  chain: string;
  resSeq: number;
  pdbResidue: string;
  dbResidue: string;
  description: string;
}

export interface PdbSite {
  name: string;
  residues: string[];
}

export interface PdbMissingResidue {
  chain: string;
  resName: string;
  resSeq: number;
}

export interface PdbLigand {
  id: string;
  chain?: string;
  name?: string;
  formula?: string;
}

export interface PdbSecondaryStructure {
  helices: Array<{id: string; chain: string; startRes: number; endRes: number; length: number}>;
  sheets: Array<{id: string; chain: string; startRes: number; endRes: number; length: number}>;
}

export interface PdbModifiedResidue {
  resName: string;
  chain: string;
  resSeq: number;
  standardRes: string;
  comment: string;
}

export interface PdbCitation {
  title?: string;
  journal?: string;
  doi?: string;
  pmid?: number;
}

export interface PdbHeaderInfo {
  pdbId?: string;
  title?: string;
  classification?: string;
  depositDate?: string;
  method?: string;
  resolution?: number;
  rFactor?: number;
  organism?: string;
  expressionSystem?: string;
  software?: string;
  softwareRemarks?: string[];
  entities?: PdbEntity[];
  mutations?: PdbMutation[];
  sites?: PdbSite[];
  missingResidues?: PdbMissingResidue[];
  citation?: PdbCitation;
  ligands?: PdbLigand[];
  secondaryStructure?: PdbSecondaryStructure;
  modifiedResidues?: PdbModifiedResidue[];
  disulfideBondCount?: number;
  atomCount?: number;
  chains?: string[];
  residueCountByChain?: {[chain: string]: number};
}

function _collectLines(lines: string[], recordName: string): string {
  const prefix = recordName.padEnd(6);
  const parts: string[] = [];
  for (const line of lines) {
    if (line.startsWith(prefix))
      parts.push(line.substring(10).trim());
  }
  return parts.join(' ').replace(/\s+/g, ' ').trim();
}

function _extractKeyValue(lines: string[], recordName: string, key: string): string | undefined {
  const fullText = _collectLines(lines, recordName);
  const regex = new RegExp(`${key}:\\s*([^;]+)`, 'i');
  const match = fullText.match(regex);
  return match ? match[1].trim() : undefined;
}

function _parseEntities(lines: string[]): PdbEntity[] | undefined {
  const fullText = _collectLines(lines, 'COMPND');
  if (!fullText) return undefined;
  const molBlocks = fullText.split(/MOL_ID:\s*/i).filter(Boolean);
  const entities: PdbEntity[] = [];
  for (const block of molBlocks) {
    const molIdMatch = block.match(/^(\d+)\s*;/);
    const molId = molIdMatch ? molIdMatch[1] : (entities.length + 1).toString();
    const moleculeMatch = block.match(/MOLECULE:\s*([^;]+)/i);
    const chainMatch = block.match(/CHAIN:\s*([^;]+)/i);
    const engineered = /ENGINEERED:\s*YES/i.test(block);
    const mutation = /MUTATION:\s*YES/i.test(block);
    entities.push({
      molId,
      molecule: moleculeMatch ? moleculeMatch[1].trim() : undefined,
      chains: chainMatch ? chainMatch[1].split(',').map((c) => c.trim()).filter(Boolean) : undefined,
      engineered: engineered || undefined,
      mutation: mutation || undefined,
    });
  }
  return entities.length > 0 ? entities : undefined;
}

function _parseCitation(lines: string[]): PdbCitation | undefined {
  const jrnlLines = lines.filter((l) => l.startsWith('JRNL  '));
  if (jrnlLines.length === 0) return undefined;
  const titleParts: string[] = [];
  const refParts: string[] = [];
  let doi: string | undefined;
  let pmid: number | undefined;
  for (const line of jrnlLines) {
    const subtype = line.substring(12, 16).trim();
    const content = line.substring(19).trim();
    if (subtype === 'TITL')
      titleParts.push(content);
    else if (subtype === 'REF')
      refParts.push(content);
    else if (subtype === 'DOI')
      doi = content;
    else if (subtype === 'PMID')
      pmid = parseInt(content, 10) || undefined;
  }
  const title = titleParts.join(' ').replace(/\s+/g, ' ').trim() || undefined;
  const journal = refParts.length > 0 ?
    refParts[0].replace(/\s{2,}/g, ' ').trim().split(/\s{2,}|V\.\s/)[0].trim() :
    undefined;
  if (!title && !journal && !doi && !pmid) return undefined;
  return {title, journal, doi, pmid};
}

function _parseLigands(lines: string[]): PdbLigand[] | undefined {
  const ligandMap: {[id: string]: PdbLigand} = {};
  for (const line of lines) {
    if (line.startsWith('HET   ')) {
      const id = line.substring(7, 10).trim();
      const chain = line.substring(12, 13).trim() || undefined;
      if (id && id !== 'HOH') {
        if (!ligandMap[id])
          ligandMap[id] = {id, chain};
      }
    }
  }
  for (const line of lines) {
    if (line.startsWith('HETNAM')) {
      const id = line.substring(11, 14).trim();
      const namePart = line.substring(15).trim();
      if (id && ligandMap[id]) {
        ligandMap[id].name = ligandMap[id].name ?
          ligandMap[id].name + ' ' + namePart : namePart;
      }
    }
  }
  for (const line of lines) {
    if (line.startsWith('FORMUL')) {
      const id = line.substring(12, 15).trim();
      const formulaPart = line.substring(18).trim();
      if (id && ligandMap[id]) {
        ligandMap[id].formula = ligandMap[id].formula ?
          ligandMap[id].formula + ' ' + formulaPart : formulaPart;
      }
    }
  }
  const ligands = Object.values(ligandMap);
  return ligands.length > 0 ? ligands : undefined;
}

function _parseMutations(lines: string[]): PdbMutation[] | undefined {
  const mutations: PdbMutation[] = [];
  for (const line of lines) {
    if (line.startsWith('SEQADV') && /ENGINEERED MUTATION/i.test(line)) {
      const pdbResidue = line.substring(12, 15).trim();
      const chain = line.substring(16, 17).trim();
      const resSeq = parseInt(line.substring(18, 22).trim(), 10);
      const dbResidue = line.substring(39, 42).trim();
      const seqNum = line.substring(43, 48).trim();
      if (chain && !isNaN(resSeq)) {
        mutations.push({
          chain, resSeq, pdbResidue, dbResidue,
          description: `${dbResidue}${seqNum}${pdbResidue}`,
        });
      }
    }
  }
  return mutations.length > 0 ? mutations : undefined;
}

function _parseSites(lines: string[]): PdbSite[] | undefined {
  const siteMap: {[name: string]: string[]} = {};
  for (const line of lines) {
    if (line.startsWith('SITE  ')) {
      const name = line.substring(11, 14).trim();
      if (!name) continue;
      if (!siteMap[name]) siteMap[name] = [];
      for (let col = 18; col < 62; col += 11) {
        const resName = line.substring(col, col + 3).trim();
        const chain = line.substring(col + 4, col + 5).trim();
        const resSeq = line.substring(col + 5, col + 9).trim();
        if (resName && chain)
          siteMap[name].push(`${resName} ${chain}${resSeq}`);
      }
    }
  }
  const sites = Object.entries(siteMap).map(([name, residues]) => ({name, residues}));
  return sites.length > 0 ? sites : undefined;
}

function _parseMissingResidues(lines: string[]): PdbMissingResidue[] | undefined {
  const missing: PdbMissingResidue[] = [];
  let inData = false;
  for (const line of lines) {
    if (!line.startsWith('REMARK 465')) continue;
    const content = line.substring(10).trim();
    if (content.includes('M RES C SSSEQI')) { inData = true; continue; }
    if (!inData || !content) continue;
    const parts = content.trim().split(/\s+/);
    if (parts.length >= 3) {
      const idx = parts.length >= 4 && /^\d+$/.test(parts[0]) ? 1 : 0;
      const resName = parts[idx];
      const chain = parts[idx + 1];
      const resSeq = parseInt(parts[idx + 2], 10);
      if (resName && chain && !isNaN(resSeq))
        missing.push({chain, resName, resSeq});
    }
  }
  return missing.length > 0 ? missing : undefined;
}

function _parseSoftware(lines: string[]): {software?: string; remarks?: string[]} {
  let software: string | undefined;
  const remarks: string[] = [];
  for (const line of lines) {
    if (line.startsWith('REMARK   3') && /PROGRAM\s+:/.test(line)) {
      const match = line.match(/PROGRAM\s*:\s*(.+)/);
      if (match) software = match[1].trim();
    }
    if (line.startsWith('REMARK  99')) {
      const content = line.substring(10).trim();
      if (content && !content.startsWith('REMARK'))
        remarks.push(content);
    }
  }
  return {software, remarks: remarks.length > 0 ? remarks : undefined};
}

function _parseSecondaryStructure(lines: string[]): PdbSecondaryStructure | undefined {
  const helices: PdbSecondaryStructure['helices'] = [];
  const sheets: PdbSecondaryStructure['sheets'] = [];
  for (const line of lines) {
    if (line.startsWith('HELIX ')) {
      const id = line.substring(11, 14).trim();
      const chain = line.substring(19, 20).trim();
      const startRes = parseInt(line.substring(21, 25).trim(), 10);
      const endRes = parseInt(line.substring(33, 37).trim(), 10);
      const length = parseInt(line.substring(71, 76).trim(), 10);
      if (chain && !isNaN(startRes) && !isNaN(endRes)) {
        helices.push({
          id, chain, startRes, endRes,
          length: !isNaN(length) ? length : endRes - startRes + 1,
        });
      }
    }
    if (line.startsWith('SHEET ')) {
      const id = line.substring(11, 14).trim();
      const chain = line.substring(21, 22).trim();
      const startRes = parseInt(line.substring(22, 26).trim(), 10);
      const endRes = parseInt(line.substring(33, 37).trim(), 10);
      if (chain && !isNaN(startRes) && !isNaN(endRes)) {
        sheets.push({
          id, chain, startRes, endRes,
          length: endRes - startRes + 1,
        });
      }
    }
  }
  if (helices.length === 0 && sheets.length === 0) return undefined;
  return {helices, sheets};
}

function _parseModifiedResidues(lines: string[]): PdbModifiedResidue[] | undefined {
  const mods: PdbModifiedResidue[] = [];
  for (const line of lines) {
    if (line.startsWith('MODRES')) {
      const resName = line.substring(12, 15).trim();
      const chain = line.substring(16, 17).trim();
      const resSeq = parseInt(line.substring(18, 22).trim(), 10);
      const standardRes = line.substring(24, 27).trim();
      const comment = line.substring(29).trim();
      if (resName && chain && !isNaN(resSeq))
        mods.push({resName, chain, resSeq, standardRes, comment});
    }
  }
  return mods.length > 0 ? mods : undefined;
}

/** Parse all available header information from a PDB file string. */
export function parsePdbHeaders(pdbText: string): PdbHeaderInfo {
  const lines = pdbText.split('\n');
  const result: PdbHeaderInfo = {};

  // HEADER: classification (cols 11-50), date (cols 51-59), PDB ID (cols 63-66)
  const headerLine = lines.find((l) => l.startsWith('HEADER'));
  if (headerLine) {
    const classification = headerLine.substring(10, 50).trim();
    if (classification) result.classification = classification;
    const dateStr = headerLine.substring(50, 59).trim();
    if (dateStr && /\d{2}-[A-Z]{3}-\d{2}/.test(dateStr)) result.depositDate = dateStr;
    const pdbId = headerLine.substring(62, 66).trim();
    if (pdbId && /^[A-Za-z0-9]{4}$/.test(pdbId)) result.pdbId = pdbId;
  }

  // TITLE
  const title = _collectLines(lines, 'TITLE ');
  if (title) result.title = title;

  // EXPDTA
  const method = _collectLines(lines, 'EXPDTA');
  if (method) result.method = method;

  // REMARK 2: resolution
  for (const line of lines) {
    if (line.startsWith('REMARK   2') && line.includes('RESOLUTION')) {
      const match = line.match(/(\d+\.\d+)\s+ANGSTROMS/);
      if (match) result.resolution = parseFloat(match[1]);
      break;
    }
  }

  // REMARK 3: R-factor
  for (const line of lines) {
    if (line.startsWith('REMARK   3') && /R VALUE\s+\(WORKING SET\)/i.test(line)) {
      const match = line.match(/:\s*(\d+\.\d+)/);
      if (match) result.rFactor = parseFloat(match[1]);
      break;
    }
  }

  // SOURCE: organism, expression system
  const organism = _extractKeyValue(lines, 'SOURCE', 'ORGANISM_SCIENTIFIC');
  if (organism) result.organism = organism;
  const sourceText = _collectLines(lines, 'SOURCE');
  if (sourceText) {
    const exprMatch = sourceText.match(/EXPRESSION_SYSTEM:\s*([^;]+)/i);
    if (exprMatch)
      result.expressionSystem = exprMatch[1].trim();
  }

  // COMPND: entities
  result.entities = _parseEntities(lines);

  // JRNL: citation
  result.citation = _parseCitation(lines);

  // MODRES: modified residues (parse before ligands to filter them out)
  result.modifiedResidues = _parseModifiedResidues(lines);
  const modResIds = new Set((result.modifiedResidues ?? []).map((m) => m.resName));

  // Ligands: HET, HETNAM, FORMUL (exclude modified residues)
  const allLigands = _parseLigands(lines);
  if (allLigands)
    result.ligands = allLigands.filter((l) => !modResIds.has(l.id));
  if (result.ligands && result.ligands.length === 0)
    result.ligands = undefined;

  // SEQADV: mutations
  result.mutations = _parseMutations(lines);

  // SITE: binding/active sites
  result.sites = _parseSites(lines);

  // REMARK 465: missing residues
  result.missingResidues = _parseMissingResidues(lines);

  // Software: REMARK 3 (refinement program) and REMARK 99 (MOE/Schrödinger)
  const sw = _parseSoftware(lines);
  if (sw.software) result.software = sw.software;
  if (sw.remarks) result.softwareRemarks = sw.remarks;

  // HELIX/SHEET: secondary structure
  result.secondaryStructure = _parseSecondaryStructure(lines);

  // SSBOND count
  const ssbondCount = lines.filter((l) => l.startsWith('SSBOND')).length;
  if (ssbondCount > 0) result.disulfideBondCount = ssbondCount;

  // Atom count, chains, and residue count per chain from ATOM records
  let atomCount = 0;
  const chainSet = new Set<string>();
  const residuesByChain: {[chain: string]: Set<string>} = {};
  for (const line of lines) {
    if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
      atomCount++;
      const chain = line.substring(21, 22).trim();
      if (chain) chainSet.add(chain);
    }
    if (line.startsWith('ATOM  ')) {
      const chain = line.substring(21, 22).trim();
      const resSeq = line.substring(22, 27).trim();
      if (chain) {
        if (!residuesByChain[chain]) residuesByChain[chain] = new Set();
        residuesByChain[chain].add(resSeq);
      }
    }
  }
  if (atomCount > 0) result.atomCount = atomCount;
  if (chainSet.size > 0) result.chains = [...chainSet].sort();
  const resCounts: {[chain: string]: number} = {};
  for (const [chain, residues] of Object.entries(residuesByChain))
    resCounts[chain] = residues.size;
  if (Object.keys(resCounts).length > 0) result.residueCountByChain = resCounts;

  return result;
}


/** {@link https://molstar.org/docs/plugin/#plugincontext-without-built-in-react-ui} */
const MolstarPluginSpec: PluginSpec = {
  // eslint-disable-next-line new-cap
  ...DefaultPluginSpec(),
  config: [
    [PluginConfig.VolumeStreaming.Enabled, false],
  ],
};


export class PdbResDataFrame extends DG.DataFrame implements PdbResDataFrameType {
  public static ColNames = {
    code: 'code',
    compId: 'compId',
    seqId: 'seqId',
    label: 'label',
    seq: 'seq',
    frame: 'frame',
  };

  public readonly code: DG.Column<string>;
  public readonly compId: DG.Column<string>;
  public readonly seqId: DG.Column<number>;
  public readonly label: DG.Column<string>;
  public readonly seq: DG.Column<string>;
  public readonly frame: DG.Column<number>;

  protected constructor(df: DG.DataFrame) {
    super(df.dart);

    this.code = this.getCol(PdbResDataFrame.ColNames.code);
    this.compId = this.getCol(PdbResDataFrame.ColNames.compId);
    this.seqId = this.getCol(PdbResDataFrame.ColNames.seqId);
    this.label = this.getCol(PdbResDataFrame.ColNames.label);
    this.seq = this.getCol(PdbResDataFrame.ColNames.seq);
    this.frame = this.getCol(PdbResDataFrame.ColNames.frame);
  }

  public static createDf(length: number,
    code: (i: number) => string, compId: (i: number) => string,
    seqId: (i: number) => number, label: (i: number) => string,
    seq: string, frame: number,
  ): PdbResDataFrame {
    const codeCol = DG.Column.string('code', length).init((i) => code(i));
    const compIdCol = DG.Column.int('compId', length).init((i) => compId(i));
    const seqIdCol = DG.Column.int('seqId', length).init((i) => seqId(i));
    const labelCol = DG.Column.string('label', length).init((i) => label(i));
    const seqCol = DG.Column.string('seq', length).init((_i) => seq);
    const frameCol = DG.Column.int('frame', length).init((_i) => frame);
    const cols = [
      codeCol,
      compIdCol,
      seqIdCol,
      labelCol,
      seqCol,
      frameCol,
    ];
    const resDf = new PdbResDataFrame(DG.DataFrame.fromColumns(cols));
    return resDf;
  }
}

type PdbHelperWindowType = Window & {
  $pdbHelper?: PdbHelper,
};
declare const window: PdbHelperWindowType;

export class PdbHelper implements IPdbHelper {
  //private stage: NGL.Stage;
  private plugin: PluginContext;

  /** Protect constructor to prevent multiple instantiation. */
  protected constructor() {
    //this.stage = new NGL.Stage();
  }

  protected async init(): Promise<void> {
    this.plugin = new PluginContext(MolstarPluginSpec);
    await this.plugin.init();
  }

  /** Parse PDB file header metadata (TITLE, SOURCE, JRNL, SITE, etc.)
   *  Complementary to Molstar's parser which handles coordinates/sequences
   *  but not header metadata. */
  static parseHeaders(pdbStr: string): PdbHeaderInfo {
    return parsePdbHeaders(pdbStr); // defined above the class
  }

  async pdbToDf(pdbStr: string, _name: string): Promise<PdbResDataFrameType> {
    //https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    // using molstar parser
    if (!pdbStr) throw new Error('Empty PDB data');

    const pdbData: StateObjectSelector = await this.plugin.builders.data.rawData({data: pdbStr});
    const trState = await this.plugin.builders.structure.parseTrajectory(pdbData, 'pdb');

    if (!trState.isOk || !trState.obj) throw new Error(`Trajectory is not Ok`);
    const tr: PluginStateObject.Molecule.Trajectory = trState.obj;

    function seqToDf(src: Sequence, seq: string, frame: number): PdbResDataFrame {
      const codeAL: ArrayLike<string> = src.code.toArray();
      const compIdAL: ArrayLike<string> = src.compId.toArray();
      const seqIdAL: ArrayLike<number> = src.seqId.toArray();
      const labelAL: ArrayLike<string> = src.label.toArray();
      const seqResDf: PdbResDataFrame = PdbResDataFrame.createDf(src.length,
        (i) => codeAL[i], (i) => compIdAL[i], (i) => seqIdAL[i], (i) => labelAL[i],
        seq, frame);
      return seqResDf;
    }

    const resDfList: DG.DataFrame[] = wu.count(0).take(tr.data.frameCount).map((frameI) => {
      const frame: Model = <Model>tr.data.getFrameAtIndex(frameI);
      const frameResDfList: DG.DataFrame[] = frame.sequence.sequences.map((seq) => {
        const seqDf: DG.DataFrame = seqToDf(seq.sequence, seq.entityId, frameI);
        return seqDf;
      });
      return frameResDfList;
    }).toArray().flat();

    const resDf: PdbResDataFrameType = PdbResDataFrame.createDf(0,
      (_i) => '', (_i) => '', (_i) => -1, (_i) => '',
      '', -1);
    for (const seqDf of resDfList)
      resDf.append(seqDf, true);

    resDf.setTag(pdbTAGS.PDB, pdbStr);
    resDf.temp.set(pdbTAGS.PDB, pdbData);

    return resDf;
  }

  parsePdbqt(pdbqtStr: string, molColName?: string): DG.DataFrame {
    const data: Pdbqt = Pdbqt.parse(pdbqtStr);
    const molColNameVal: string = molColName ?? IMPORT[Molecule3DUnits.pdbqt].molColName;
    const resDf = data.toDataFrame(molColNameVal);
    return resDf;
  }

  async molToPdb(mol: string): Promise<string> {
    // const val = await ngl.autoLoad(mol, {ext: 'sdf'});
    // const resPdb = (new ngl.PdbWriter(val)).getString();
    // return resPdb;
    const resName: string = 'UNK';
    const chain: string = '';
    const resNum: number = 0;

    const molH = MolfileHandler.getInstance(mol);
    const lineList: LineBase[] = new Array<LineBase>(molH.atomCount + 1);
    for (let atomI = 0; atomI < molH.atomCount; ++atomI) {
      const atomType = molH.atomTypes[atomI];
      const atomX: number = molH.x[atomI];
      const atomY: number = molH.y[atomI];
      const atomZ: number = molH.z[atomI];
      // @formatter:off
      lineList[atomI] = new PdbAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('HETATM'),
        (atomI + 1), atomType, '', '', resName, chain, resNum, ''),
      atomX, atomY, atomZ, 0, 0), '', atomType, '');
      // @formatter:on
    }
    lineList[molH.atomCount] = new PdbAtomTer(new AtomBase(new LineBase('TER'),
      -1, '', '', '', '', '', -1, ''));

    const resPdb = lineList.map((l) => l.toStr()).join('\n');
    return resPdb;
  }

  async pdbqtToMol(srcPdbqt: string): Promise<string> {
    const srcBlob = new Blob([srcPdbqt]);
    const valS: ngl.Structure = await ngl.autoLoad(srcBlob, {ext: 'pdbqt'});

    const res = (new ngl.SdfWriter(valS)).getData();
    return res;
  }

  // -- Instance singleton --

  public static async getInstance(): Promise<IPdbHelper> {
    let res: PdbHelper = window.$pdbHelper!;
    if (!res) {
      window.$pdbHelper = res = new PdbHelper();
      await res.init();
    }
    return res;
  }
}
