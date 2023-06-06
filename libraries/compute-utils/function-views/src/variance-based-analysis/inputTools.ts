// inputTools.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getGeneratedColumnsData} from './inputsGeneration';

export type VariedNumericalInputInfo = {
  prop: DG.Property,
  min: number,
  max: number
};

export type FixedInputItem = {
  name: string,
  value: any
};

function getNumericalInputCol(numInputInfo: VariedNumericalInputInfo,
  randomData: Float32Array): DG.Column {
  const length = randomData.length;

  const column = DG.Column.fromType(numInputInfo.prop.propertyType as unknown as DG.COLUMN_TYPE,
    numInputInfo.prop.caption ?? numInputInfo.prop.name, length);

  const columnSource = column.getRawData();

  const step = numInputInfo.max - numInputInfo.min;

  for (let i = 0; i < length; ++i)
    columnSource[i] = numInputInfo.min + step * randomData[i];

  return column;
}

export function getVariedNumericalInputColumns(samplesCount: number,
  inputsInfo: VariedNumericalInputInfo[]): DG.Column[] {
  const dimension = inputsInfo.length;

  const randData = getGeneratedColumnsData(samplesCount, dimension);

  return [...Array(dimension).keys()].map((i) => getNumericalInputCol(inputsInfo[i], randData[i]));
}
