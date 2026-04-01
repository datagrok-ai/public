/** Pure PDB file header parser. All fields are optional to handle
 *  files from RCSB, MOE, Schrödinger, and other sources. */

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

/** Collect continuation lines for a given record type (e.g. TITLE, EXPDTA).
 *  PDB continuation lines have the same record name with a serial number in cols 9-10. */
function collectLines(lines: string[], recordName: string): string {
  const prefix = recordName.padEnd(6);
  const parts: string[] = [];
  for (const line of lines) {
    if (line.startsWith(prefix))
      parts.push(line.substring(10).trim());
  }
  return parts.join(' ').replace(/\s+/g, ' ').trim();
}

/** Extract a value from SOURCE/COMPND style key-value lines.
 *  These records use "KEY: value;" format across continuation lines. */
function extractKeyValue(lines: string[], recordName: string, key: string): string | undefined {
  const fullText = collectLines(lines, recordName);
  // Match KEY: value; or KEY: value (at end of string)
  const regex = new RegExp(`${key}:\\s*([^;]+)`, 'i');
  const match = fullText.match(regex);
  return match ? match[1].trim() : undefined;
}

/** Parse COMPND records into entities grouped by MOL_ID. */
function parseEntities(lines: string[]): PdbEntity[] | undefined {
  const fullText = collectLines(lines, 'COMPND');
  if (!fullText) return undefined;

  // Split by MOL_ID
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

/** Parse JRNL records into citation info. */
function parseCitation(lines: string[]): PdbCitation | undefined {
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

/** Parse HET/HETNAM/FORMUL records into ligand info. */
function parseLigands(lines: string[]): PdbLigand[] | undefined {
  const ligandMap: {[id: string]: PdbLigand} = {};

  // HET records: residue name, chain, number of HETATM atoms
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

  // HETNAM records: compound name
  for (const line of lines) {
    if (line.startsWith('HETNAM')) {
      const id = line.substring(11, 14).trim();
      const namePart = line.substring(15).trim();
      if (id && ligandMap[id]) {
        ligandMap[id].name = ligandMap[id].name ?
          ligandMap[id].name + ' ' + namePart :
          namePart;
      }
    }
  }

  // FORMUL records: chemical formula
  for (const line of lines) {
    if (line.startsWith('FORMUL')) {
      const id = line.substring(12, 15).trim();
      const formulaPart = line.substring(18).trim();
      if (id && ligandMap[id]) {
        ligandMap[id].formula = ligandMap[id].formula ?
          ligandMap[id].formula + ' ' + formulaPart :
          formulaPart;
      }
    }
  }

  const ligands = Object.values(ligandMap);
  return ligands.length > 0 ? ligands : undefined;
}

/** Parse SEQADV records for engineered mutations. */
function parseMutations(lines: string[]): PdbMutation[] | undefined {
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
          chain,
          resSeq,
          pdbResidue,
          dbResidue,
          description: `${dbResidue}${seqNum}${pdbResidue}`,
        });
      }
    }
  }
  return mutations.length > 0 ? mutations : undefined;
}

/** Parse SITE records into binding/active site definitions. */
function parseSites(lines: string[]): PdbSite[] | undefined {
  const siteMap: {[name: string]: string[]} = {};
  for (const line of lines) {
    if (line.startsWith('SITE  ')) {
      const name = line.substring(11, 14).trim();
      if (!name) continue;
      if (!siteMap[name]) siteMap[name] = [];
      // Residues are in groups of 4, starting at column 18
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

/** Parse REMARK 465 for missing residues. */
function parseMissingResidues(lines: string[]): PdbMissingResidue[] | undefined {
  const missing: PdbMissingResidue[] = [];
  let inData = false;
  for (const line of lines) {
    if (!line.startsWith('REMARK 465')) continue;
    const content = line.substring(10).trim();
    // Skip header lines; data starts after "M RES C SSSEQI" line
    if (content.includes('M RES C SSSEQI')) {
      inData = true;
      continue;
    }
    if (!inData || !content) continue;
    // Format: [model] resName chain resSeq [iCode]
    const parts = content.trim().split(/\s+/);
    if (parts.length >= 3) {
      const idx = parts.length >= 4 && /^\d+$/.test(parts[0]) ? 1 : 0; // skip model number
      const resName = parts[idx];
      const chain = parts[idx + 1];
      const resSeq = parseInt(parts[idx + 2], 10);
      if (resName && chain && !isNaN(resSeq))
        missing.push({chain, resName, resSeq});
    }
  }
  return missing.length > 0 ? missing : undefined;
}

/** Extract software info from REMARK 3 (refinement) and REMARK 99 (MOE/Schrödinger). */
function parseSoftware(lines: string[]): {software?: string; remarks?: string[]} {
  let software: string | undefined;
  const remarks: string[] = [];

  for (const line of lines) {
    // REMARK 3: refinement program
    if (line.startsWith('REMARK   3') && /PROGRAM\s+:/.test(line)) {
      const match = line.match(/PROGRAM\s*:\s*(.+)/);
      if (match) software = match[1].trim();
    }
    // REMARK 99: MOE, Schrödinger, etc.
    if (line.startsWith('REMARK  99')) {
      const content = line.substring(10).trim();
      if (content && !content.startsWith('REMARK'))
        remarks.push(content);
    }
  }

  return {software, remarks: remarks.length > 0 ? remarks : undefined};
}

/** Parse HELIX and SHEET records for secondary structure. */
function parseSecondaryStructure(lines: string[]): PdbSecondaryStructure | undefined {
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

/** Parse MODRES records for modified residues. */
function parseModifiedResidues(lines: string[]): PdbModifiedResidue[] | undefined {
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
  const title = collectLines(lines, 'TITLE ');
  if (title) result.title = title;

  // EXPDTA
  const method = collectLines(lines, 'EXPDTA');
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
  const organism = extractKeyValue(lines, 'SOURCE', 'ORGANISM_SCIENTIFIC');
  if (organism) result.organism = organism;

  const exprSystem = extractKeyValue(lines, 'SOURCE', 'EXPRESSION_SYSTEM(?!_)');
  // More precise: get the value directly
  const sourceText = collectLines(lines, 'SOURCE');
  if (sourceText) {
    const exprMatch = sourceText.match(/EXPRESSION_SYSTEM:\s*([^;]+)/i);
    if (exprMatch) {
      // Make sure we don't capture EXPRESSION_SYSTEM_TAXID etc.
      result.expressionSystem = exprMatch[1].trim();
    }
  }

  // COMPND: entities
  result.entities = parseEntities(lines);

  // JRNL: citation
  result.citation = parseCitation(lines);

  // Ligands: HET, HETNAM, FORMUL
  result.ligands = parseLigands(lines);

  // SEQADV: mutations
  result.mutations = parseMutations(lines);

  // SITE: binding/active sites
  result.sites = parseSites(lines);

  // REMARK 465: missing residues
  result.missingResidues = parseMissingResidues(lines);

  // Software: REMARK 3 (refinement program) and REMARK 99 (MOE/Schrödinger)
  const sw = parseSoftware(lines);
  if (sw.software) result.software = sw.software;
  if (sw.remarks) result.softwareRemarks = sw.remarks;

  // HELIX/SHEET: secondary structure
  result.secondaryStructure = parseSecondaryStructure(lines);

  // MODRES: modified residues
  result.modifiedResidues = parseModifiedResidues(lines);

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
    // Count unique residues from ATOM records only (not HETATM)
    if (line.startsWith('ATOM  ')) {
      const chain = line.substring(21, 22).trim();
      const resSeq = line.substring(22, 27).trim(); // includes insertion code
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
