// outputTools.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const UNSUPPORTED_TYPE_ERROR_MSG = 'Unsupported output type.';

type OutputInfo = {
  name: string,  
  type: DG.TYPE
};

function getOutputsInfo(func: DG.Func): Array<OutputInfo> {
  const ouputInfo = Array<OutputInfo>(0);

  for (const item of func.outputs)
    ouputInfo.push({name: item.name, type: item.propertyType} as OutputInfo);

  return ouputInfo;
}

function getDataFrameItemsInfo(funcCall: DG.FuncCall, outputName: string): Array<OutputInfo> {
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

export function getEachOutputItemInfo(funcCall: DG.FuncCall): Array<OutputInfo> {
  const outputInfoBasic = getOutputsInfo(funcCall.func);

  const ouputInfoExtended = Array<OutputInfo>(0);

  for (const item of outputInfoBasic) {
    console.log(item);

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

/*
const INDEX_OF_NAME = 0;
const INDEX_OF_VALUE = 1;

export function getOutputTypes(func: DG.Func): Array<DG.TYPE> {
  const outputTypes: Array<DG.TYPE> = [];

  for (const item of func.outputs)
    outputTypes.push(item.propertyType);

  return outputTypes;
}



function getOutputItemsProperyInfo(outputType: DG.TYPE, output: any): Array<PropertyInfo> {
  const ouputInfo = Array<PropertyInfo>(0);

  console.log(outputType);

  switch (outputType) {
    case DG.TYPE.BIG_INT:
    case DG.TYPE.INT:
    case DG.TYPE.FLOAT:
        ouputInfo.push({name: output[INDEX_OF_NAME], root: null} as PropertyInfo);
        break;
    case DG.TYPE.DATA_FRAME:
        ouputInfo.push.apply(output[INDEX_OF_VALUE]);
        break;
    default:
        break;
  }

  return ouputInfo;
}

export function getOutputNames(funcCall: DG.FuncCall): Array<PropertyInfo> {
  const outputTypes = getOutputTypes(funcCall.func);

  const ouputInfo = Array<PropertyInfo>(0);
  
  let index: number = 0;

  for (const output of [...funcCall.outputs]) {
    console.log(getOutputItemsProperyInfo(outputTypes[index], output));
    ++index;
  }  

  return ouputInfo;
} */