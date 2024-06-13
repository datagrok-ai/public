// Regression tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_fitLinearRegressionParams} from '../wasm/EDAAPI';

/** Computes coefficients of linear regression */
export async function computeLinRegressionCoefs(features: DG.ColumnList, targets: DG.Column): Promise<void> {
  const featuresCount = features.length;
  const paramsCount = featuresCount + 1;

  const params = _fitLinearRegressionParams(features, targets, paramsCount) as DG.Column;

  console.log(params);
  console.log(params.getRawData());
}
