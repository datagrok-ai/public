// outputTools.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const UNSUPPORTED_TYPE_ERROR_MSG = 'Unsupported output type.';

type OutputInfo = {
  name: string,  
  type: DG.TYPE
};

function getOutputsInfo(func: DG.Func): OutputInfo[] {
  return func.outputs.map(item => ({name: item.name, type: item.propertyType}));
}

function getDataFrameItemsInfo(funcCall: DG.FuncCall, outputName: string): OutputInfo[] {
  const ouputInfo = Array<OutputInfo>(0);

  const df = funcCall.outputs[outputName] as DG.DataFrame;
  
  if (df.rowCount === 1)
    for (const col of df.columns)     
      ouputInfo.push({name: `${col.name} (${outputName})`, type: col.type} as OutputInfo);
  else
    for (const col of df.columns)
      for (let i = 1; i <= df.rowCount; ++i)
        ouputInfo.push({name: `${col.name} (${outputName}, ${i})`, type: col.type} as OutputInfo);
  
  return ouputInfo;
}

export function getEachOutputItemInfo(funcCall: DG.FuncCall): OutputInfo[] {
  const outputInfoBasic = getOutputsInfo(funcCall.func);

  const ouputInfoExtended = Array<OutputInfo>(0);

  for (const item of outputInfoBasic) {
    switch (item.type) {
      case DG.TYPE.BIG_INT:
      case DG.TYPE.INT:
      case DG.TYPE.FLOAT:
        ouputInfoExtended.push(item);
        break;
      case DG.TYPE.DATA_FRAME:
        Array.prototype.push.apply(ouputInfoExtended, getDataFrameItemsInfo(funcCall, item.name));
        break;
      default:
        throw new Error(UNSUPPORTED_TYPE_ERROR_MSG);
    }
  }
  
  return ouputInfoExtended;
}

export function getOutputColumns(funcCalls: DG.FuncCall[]): DG.Column[] {
  const columns: DG.Column[] = [];

  const outputsInfoBasic = getOutputsInfo(funcCalls[0].func);
  const outputsInfoExtended = getEachOutputItemInfo(funcCalls[0]);

  const rowCount = funcCalls.length;

  for (const item of outputsInfoExtended)
    columns.push(DG.Column.fromType(item.type as unknown as DG.COLUMN_TYPE, item.name, rowCount));

  for (let row = 0; row < rowCount; ++row) {
    const funcCall = funcCalls[row];

    let col = 0;

    for (const item of outputsInfoBasic)
      switch (item.type) {
        case DG.TYPE.BIG_INT:
        case DG.TYPE.INT:
        case DG.TYPE.FLOAT:
          columns[col].set(row, funcCall.outputs[item.name]);
          ++col;
          break;
        case DG.TYPE.DATA_FRAME:
          for (const column of funcCall.outputs[item.name].columns)
            for (const element of column.getRawData()) {
              columns[col].set(row, element);
              ++col;
            }
          break;
        default:
          throw new Error(UNSUPPORTED_TYPE_ERROR_MSG);
    }
  }

  return columns;
}