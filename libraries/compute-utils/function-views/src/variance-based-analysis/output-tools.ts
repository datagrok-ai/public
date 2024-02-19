// output-tools.ts

// Tools for operating outputs for variance-based sensitivity analysis (VSA).

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const UNSUPPORTED_TYPE_ERROR_MSG = 'Unsupported output type.';

// Output item specification
type OutputInfo = {
  name: string,
  type: DG.TYPE
};

// Returns an array that contains specification of each output item
function getOutputsInfo(func: DG.Func): OutputInfo[] {
  return func.outputs.map((item) => ({name: item.name, type: item.propertyType}));
}

// Returns specificatios corresponding to each dataframe element
function getDataFrameItemsInfo(funcCall: DG.FuncCall, outputName: string): OutputInfo[] {
  const ouputInfo = Array<OutputInfo>(0);

  const df = funcCall.outputs[outputName] as DG.DataFrame;

  /* Column-by-column scan is implemented.
     Name of each element is construced as follows:
        - [column name] ([dataframe name]) in the case of one row;
        - [column name] ([dataframe name], [no. of row]) in the case of one row. */
  if (df.rowCount === 1) {
    for (const col of df.columns)
      ouputInfo.push({name: `${col.name} (${outputName})`, type: col.type} as OutputInfo);
  } else {
    for (const col of df.columns) {
      for (let i = 1; i <= df.rowCount; ++i)
        ouputInfo.push({name: `${col.name} (${outputName}, ${i})`, type: col.type} as OutputInfo);
    }
  }

  return ouputInfo;
}

// Returns outputs specifications (dataframe elements are presented separately)
export function getEachOutputItemInfo(funcCall: DG.FuncCall): OutputInfo[] {
  // get array of outputs specifications
  const outputInfoBasic = getOutputsInfo(funcCall.func);

  const ouputInfoExtended = Array<OutputInfo>(0);

  /* create extended version of outputs specifications array:
     each dataframe element is presented separately (for the purposes of VBA) */
  for (const item of outputInfoBasic) {
    switch (item.type) {
    case DG.TYPE.BIG_INT:
    case DG.TYPE.INT:
    case DG.TYPE.FLOAT:
      ouputInfoExtended.push(item);
      break;
    case DG.TYPE.DATA_FRAME:
      //Array.prototype.push.apply(ouputInfoExtended, getDataFrameItemsInfo(funcCall, item.name));
      break;
    default:
      throw new Error(UNSUPPORTED_TYPE_ERROR_MSG);
    }
  }

  return ouputInfoExtended;
}

// Returns array of columns with results of the function calls
export function getOutputColumns(funcCalls: DG.FuncCall[]): DG.Column[] {
  const columns: DG.Column[] = [];

  const outputsInfoBasic = getOutputsInfo(funcCalls[0].func);
  const outputsInfoExtended = getEachOutputItemInfo(funcCalls[0]);

  const rowCount = funcCalls.length;

  // create column for each output
  for (const item of outputsInfoExtended)
    columns.push(DG.Column.fromType(item.type as unknown as DG.COLUMN_TYPE, item.name, rowCount));

  // fill columns with outputs (row-by-row)
  for (let row = 0; row < rowCount; ++row) {
    const funcCall = funcCalls[row];

    let col = 0;

    for (const item of outputsInfoBasic) {
      switch (item.type) {
      case DG.TYPE.BIG_INT:
      case DG.TYPE.INT:
      case DG.TYPE.FLOAT:
        columns[col].set(row, funcCall.outputs[item.name]);
        ++col;
        break;
      case DG.TYPE.DATA_FRAME:
        // each element of dataframe is put to the corresponding column (its own)
        /*for (const column of funcCall.outputs[item.name].columns) {
          for (const element of column.getRawData()) {
            columns[col].set(row, element);
            ++col;
          }
        }*/
        break;
      default:
        throw new Error(UNSUPPORTED_TYPE_ERROR_MSG);
      }
    }
  }

  return columns;
}
