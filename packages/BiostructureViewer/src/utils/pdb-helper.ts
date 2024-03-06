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
