// input-tools.ts

// Tools for operating inputs for variance-based sensitivity analysis (VSA).

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getGeneratedColumnsData} from './inputs-generation';

// Varied numerical input specification
export type VariedNumericalInputInfo = {
  prop: DG.Property,
  min: number,
  max: number
};

// Fixed input specification
export type FixedInputItem = {
  name: string,
  value: any
};

// Returns column of inputs that are used in VSA
function getNumericalInputCol(numInputInfo: VariedNumericalInputInfo,
  randomData: Float32Array): DG.Column {
  const length = randomData.length;

  const column = DG.Column.fromType(numInputInfo.prop.propertyType as unknown as DG.COLUMN_TYPE,
    numInputInfo.prop.caption ?? numInputInfo.prop.name, length);

  const columnSource = column.getRawData();

  /* REMARK.
     Elements of the array randomData (uniform distribution on [0, 1]) are scaled
     to the range [min, max], where min & max are defined by numInputInfo. */

  const step = numInputInfo.max - numInputInfo.min;

  for (let i = 0; i < length; ++i)
    columnSource[i] = numInputInfo.min + step * randomData[i];

  return column;
}

// Returns columns of varied numerical inputs that are used in VSA
export function getVariedNumericalInputColumnsForSobolAnalysis(samplesCount: number,
  inputsInfo: VariedNumericalInputInfo[]): DG.Column[] {
  const dimension = inputsInfo.length;

  // Get random data with the uniform distribution on [0, 1]
  const randData = getGeneratedColumnsData(samplesCount, dimension);

  // Create & return an array of columns of varied numerical inputs
  return [...Array(dimension).keys()].map((i) => getNumericalInputCol(inputsInfo[i], randData[i]));
}

// Returns columns of varied numerical inputs that are used in RSA
export function getVariedNumericalInputColumnsForRandomAnalysis(samplesCount: number,
  inputsInfo: VariedNumericalInputInfo[]): DG.Column[] {
  const columns = [] as DG.Column[];

  for (const item of inputsInfo) {
    const column = DG.Column.fromType(item.prop.propertyType as unknown as DG.COLUMN_TYPE,
      item.prop.caption ?? item.prop.name, samplesCount);

    const columnSource = column.getRawData();
    const step = item.max - item.min;

    for (let i = 0; i < samplesCount; ++i)
      columnSource[i] = item.min + step * Math.random();

    columns.push(column);
  }

  return columns;
}
