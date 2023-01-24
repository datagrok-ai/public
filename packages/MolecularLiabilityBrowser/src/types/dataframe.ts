import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export const EMPTY_ANTIGEN_NAME = '[none]';

export class AntigenDataFrame extends DG.DataFrame {
  public readonly idCol: DG.Column;
  public readonly antigenCol: DG.Column;
  public readonly antigenGsCol: DG.Column;
  public readonly clonesCol: DG.Column;
  public readonly vrsCol: DG.Column;

  protected constructor(dart: any) {
    super(dart);

    this.idCol = this.getCol('id');
    this.antigenCol = this.getCol('antigen');
    this.antigenGsCol = this.getCol('antigen_gene_symbol');
    this.clonesCol = this.getCol('clones');
    this.vrsCol = this.getCol('vrs');
  }

  public static wrap(df: DG.DataFrame): AntigenDataFrame {
    return new AntigenDataFrame(df.dart);
  }
}

/**  */
export class MlbDataFrame extends DG.DataFrame {
  public readonly antigen: string;
  public readonly vIdCol: DG.Column;

  protected constructor(dart: any, antigen: string) {
    super(dart);

    this.antigen = antigen;
    this.vIdCol = this.getCol('v_id');

    this.setTag('antigen', this.antigen);
  }

  public static wrap(df: DG.DataFrame, antigen: string): MlbDataFrame {
    return new MlbDataFrame(df.dart, antigen);
  }

  private static _empty: MlbDataFrame;

  public static get Empty(): MlbDataFrame {
    /*
    'v id', 'gdb id mappings',
    'cdr length', 'surface cdr hydrophobicity',
    'positive cdr charge', 'negative cdr charge', 'SFvCSP',
    'Heavy chain sequence', 'Light chain sequence', 'clones'
     */
    if (!MlbDataFrame._empty) {
      MlbDataFrame._empty = MlbDataFrame.wrap(DG.DataFrame.fromColumns([
        DG.Column.fromStrings('v_id', []),
        DG.Column.fromStrings('gdb_id_mappings', []),

        DG.Column.fromFloat32Array('cdr_length', new Float32Array(0)),
        DG.Column.fromFloat32Array('surface_cdr_hydrophobicity', new Float32Array(0)),
        DG.Column.fromFloat32Array('positive_cdr_charge', new Float32Array(0)),
        DG.Column.fromFloat32Array('negative_cdr_charge', new Float32Array(0)),
        DG.Column.fromFloat32Array('SFvCSP', new Float32Array(0)),

        DG.Column.fromStrings('Heavy chain sequence', []),
        DG.Column.fromStrings('Light chain sequence', []),
        DG.Column.fromStrings('clones', []), // TODO: Check index for request
        // DG.Column.from
      ]), EMPTY_ANTIGEN_NAME /* MlbDataFrame.Empty */);
    }
    return MlbDataFrame._empty;
  }
}

/** */
export class TreeDataFrame extends DG.DataFrame {
  public readonly antigen: string;
  public readonly treeCol: DG.Column;
  public readonly cloneCol: DG.Column;

  protected constructor(dart: any, antigen: string) {
    super(dart);

    this.antigen = antigen;
    this.treeCol = this.getCol('TREE');
    this.cloneCol = this.getCol('CLONE');

    this.setTag('antigen', this.antigen);
  }

  public static wrap(df: DG.DataFrame, antigen: string): TreeDataFrame {
    return new TreeDataFrame(df.dart, antigen);
  }

  private static _empty: TreeDataFrame;

  public static get Empty(): TreeDataFrame {
    /*
      'tree_id', 'antigen_id', 'antigen', 'CLONE', 'NSEQ',
      'NSITE', 'TREE_LENGTH', 'LHOOD', 'KAPPA_MLE', 'OMEGA_FWR_MLE', 'OMEGA_CDR_MLE',
      'WRC_2_MLE', 'GYW_0_MLE', 'WA_1_MLE', 'TW_0_MLE', 'SYC_2_MLE', 'GRS_0_MLE',
      'TREE'*/
    if (!TreeDataFrame._empty) {
      TreeDataFrame._empty = TreeDataFrame.wrap(DG.DataFrame.fromColumns([
        DG.Column.fromInt32Array('tree_id', new Int32Array(0)),
        DG.Column.fromInt32Array('antigen_id', new Int32Array(0)),
        DG.Column.fromStrings('antigen', []),
        DG.Column.fromStrings('CLONE', []),
        DG.Column.fromInt32Array('NSEQ', new Int32Array(0)),

        DG.Column.fromFloat32Array('NSITE', new Float32Array(0)),
        DG.Column.fromFloat32Array('TREE_LENGTH', new Float32Array(0)),
        DG.Column.fromFloat32Array('LHOOD', new Float32Array(0)),
        DG.Column.fromFloat32Array('KAPPA_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('OMEGA_FWR_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('OMEGA_CDR_MLE', new Float32Array(0)),

        DG.Column.fromFloat32Array('WRC_2_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('GYW_0_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('WA_1_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('TW_0_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('SYC_2_MLE', new Float32Array(0)),
        DG.Column.fromFloat32Array('GRS_0_MLE', new Float32Array(0)),

        DG.Column.fromStrings('TREE', []),
      ]), EMPTY_ANTIGEN_NAME /* AntigenTreeDataFrame.Empty */);
    }
    return TreeDataFrame._empty;
  }
}
