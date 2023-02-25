import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {IPdbHelper, PdbResDataFrameType} from '@datagrok-libraries/bio/src/pdb/pdb-helper';

import * as NGL from 'NGL';
import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {DefaultPluginSpec, PluginSpec} from 'molstar/lib/mol-plugin/spec';
import {PluginConfig} from 'molstar/lib/mol-plugin/config';
import {StateObject, StateObjectSelector} from 'molstar/lib/mol-state';
import {Model, Trajectory} from 'molstar/lib/mol-model/structure';
import {PluginStateObject} from 'molstar/lib/mol-plugin-state/objects';
import {Sequence} from 'molstar/lib/mol-model/sequence';

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
    frame: 'frame'
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
    seq: string, frame: number
  ): PdbResDataFrame {
    const codeCol = DG.Column.string('code', length).init((i) => code(i));
    const compIdCol = DG.Column.int('compId', length).init((i) => compId(i));
    const seqIdCol = DG.Column.int('seqId', length).init((i) => seqId(i));
    const labelCol = DG.Column.string('label', length).init((i) => label(i));
    const seqCol = DG.Column.string('seq', length).init((i) => seq);
    const frameCol = DG.Column.int('frame', length).init((i) => frame);
    const cols = [
      codeCol,
      compIdCol,
      seqIdCol,
      labelCol,
      seqCol,
      frameCol
    ];
    const resDf = new PdbResDataFrame(DG.DataFrame.fromColumns(cols));
    return resDf;
  }
}

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

  async pdbToDf(pdbStr: string, name: string): Promise<PdbResDataFrameType> {
    //https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    // using molstar parser
    const pdbData: StateObjectSelector = await this.plugin.builders.data.rawData({data: pdbStr});
    const trState = await this.plugin.builders.structure.parseTrajectory(pdbData, 'pdb');

    if (!trState.isOk || !trState.obj) throw new Error(`Trajectory is not Ok`);
    const tr: PluginStateObject.Molecule.Trajectory = trState.obj;
    const trLabel: string = tr.label;

    function seqToDf(src: Sequence, seq: string, frame: number): PdbResDataFrame {
      const length: number = src.length;
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
      (i) => '', (i) => '', (i) => -1, (i) => '',
      '', -1);
    for (const seqDf of resDfList)
      resDf.append(seqDf, true);

    resDf.setTag(pdbTAGS.PDB, pdbStr);
    resDf.temp.set(pdbTAGS.PDB, pdbData);

    return resDf;
  }

  // -- Instance singleton --
  private static _instance: PdbHelper | null;

  public static async getInstance(): Promise<IPdbHelper> {
    if (!PdbHelper._instance) {
      PdbHelper._instance = new PdbHelper();
      await PdbHelper._instance.init();
    }
    return PdbHelper._instance;
  }
}
