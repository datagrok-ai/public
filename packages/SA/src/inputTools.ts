// inputTools.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getGeneratedColumnsData} from './inputsGeneration';

export type VariedNumericalInputInfo = {
  name: string,
  caption: string | undefined,
  type: DG.COLUMN_TYPE,
  min: number,
  max: number,
  column: DG.Column | undefined
};

export type FixedInputItem = {
  name: string,
  value: any
};

function createNumericalInputCol(numInputInfo: VariedNumericalInputInfo,
  randomData: Float32Array): DG.Column {  

  const length = randomData.length;  
  
  const column = DG.Column.fromType( numInputInfo.type, numInputInfo.caption ?? numInputInfo.name, length);
  const columnSource = column.getRawData();

  const step = numInputInfo.max - numInputInfo.min;

  for (let i = 0; i < length; ++i)
    columnSource[i] = numInputInfo.min + step * randomData[i];

  numInputInfo.column = column;
  
  return column;
}

export function createVariedNumericalInputColumns(samplesCount: number, 
  inputsInfo: Array<VariedNumericalInputInfo>): void { 

  const dimension = inputsInfo.length;

  const randData = getGeneratedColumnsData(samplesCount, dimension);

  for (let i = 0; i < dimension; ++i)
    createNumericalInputCol(inputsInfo[i], randData[i]);
}