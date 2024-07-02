// Predicitve tools based on the PLS method

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Acceptable column types */
const NUMERIC_TYPES = [
  DG.COLUMN_TYPE.INT,
  DG.COLUMN_TYPE.FLOAT,
  DG.COLUMN_TYPE.BIG_INT,
  DG.COLUMN_TYPE.QNUM,
];

/** */
type PLSmodelSpecification = {
  featureNames: string[],
  regrCoefs: Float32Array,
  loadings: Float32Array[],
}

/** */
export class PLSmodel {
  private model: PLSmodelSpecification | null = null;

  constructor(packedModel?: Uint8Array) {
    if (packedModel) {}
  }

  public async fit(features: DG.ColumnList, target: DG.Column, components: number) {}

  public toBytes(): Uint8Array {
    if (this.model === null)
      throw new Error('Packing failed: model is not trained');

    return new Uint8Array();
  }

  public predict(features: DG.ColumnList): DG.Column {
    return DG.Column.fromFloat32Array(
      'prediction',
      new Float32Array(features.byIndex(0).length),
      features.byIndex(0).length,
    );
  }

  static isApplicable(features: DG.ColumnList): boolean {
    for (const col of features) {
      if (!NUMERIC_TYPES.includes(col.type as DG.COLUMN_TYPE))
        return false;
    }

    return true;
  }
};
