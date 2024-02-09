//

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {checkSize} from './utils';

const UNSUPPORTED_TYPE_ERROR_MSG = 'Unsupported output type.';
const FUNCCALL_ARRAY_IS_EMPTY_MSG = 'Funccall array is empty.';
const INCORRECT_ROW_MSG = 'Incorrect row number.';

type CellCoordinates = {
  idx: number,
  columnName: string,
};

export type OutputInfo = {
  prop: DG.Property,
  elements: CellCoordinates [],
  row: number,
};

export type SensitivityAnalysisResult ={
  funcEvalResults: DG.DataFrame,
  funcCalls: DG.FuncCall[],
};

function getCellCoordinatesFromRow(df: DG.DataFrame, row: number): CellCoordinates[] {
  const rowCount = df.rowCount;
  let idx: number = 0;

  if ((row >= 1) && (row <= rowCount))
    idx = row - 1;
  else if (row == -1)
    idx = rowCount - 1;
  else
    throw new Error(INCORRECT_ROW_MSG);

  return df.columns.names().map((colName) => ({idx: idx, columnName: colName}));
}

function getExtendedOutputsSpecification(funcCall: DG.FuncCall, outputsSpecification: OutputInfo[]): OutputInfo[] {
  const extendedOutputsSpecification = [] as OutputInfo[];

  for (const item of outputsSpecification) {
    if (item.prop.propertyType == DG.TYPE.DATA_FRAME) {
      extendedOutputsSpecification.push({
        prop: item.prop,
        elements: getCellCoordinatesFromRow(funcCall.outputs[item.prop.name], item.row),
        row: item.row,
      });
    } else
      extendedOutputsSpecification.push(item);
  }

  return extendedOutputsSpecification;
}

function getEmptyOutputTable(funcCall: DG.FuncCall, outputsSpecification: OutputInfo[], rowCount: number): DG.DataFrame {
  const columns = [] as DG.Column[];

  for (const item of outputsSpecification) {
    const prop = item.prop;
    const name = prop.name;

    switch (prop.propertyType) {
    case DG.TYPE.BIG_INT:
    case DG.TYPE.INT:
    case DG.TYPE.FLOAT:
      columns.push(DG.Column.fromType(prop.propertyType, name, rowCount));
      break;
    case DG.TYPE.DATA_FRAME:
      const df = funcCall.outputs[name] as DG.DataFrame;

      if (item.elements) {
        for (const elem of item.elements) {
          columns.push(DG.Column.fromType(
              df.col(elem.columnName)?.type as DG.COLUMN_TYPE,
              elem.columnName,
              //`${elem.columnName} `, // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13415
              //`${name}(${elem.columnName}, ${elem.idx})`,
              rowCount,
          ));
        }
      }
      break;
    default:
      throw new Error(UNSUPPORTED_TYPE_ERROR_MSG);
    }
  }
  return DG.DataFrame.fromColumns(columns);
}

export function getOutput(funcCalls: DG.FuncCall[], outputsSpecification: OutputInfo[]): DG.DataFrame {
  if (funcCalls.length < 1)
    throw new Error(FUNCCALL_ARRAY_IS_EMPTY_MSG);

  const rowCount = funcCalls.length;

  const extendedOutputsSpecification = getExtendedOutputsSpecification(funcCalls[0], outputsSpecification);

  const table = getEmptyOutputTable(funcCalls[0], extendedOutputsSpecification, rowCount);

  for (let row = 0; row < rowCount; ++row) {
    for (const item of extendedOutputsSpecification) {
      const prop = item.prop;
      const name = prop.name;
      const funcCall = funcCalls[row];

      switch (prop.propertyType) {
      case DG.TYPE.BIG_INT:
      case DG.TYPE.INT:
      case DG.TYPE.FLOAT:
        table.col(name)?.set(row, funcCall.outputs[name]);
        break;

      case DG.TYPE.DATA_FRAME:
        const df = funcCall.outputs[name] as DG.DataFrame;

        if (item.elements) {
          for (const elem of item.elements) {
            table.col(elem.columnName)?.set(
              row,
              df.col(elem.columnName)?.get(elem.idx),
            );
          }
        }

        break;
      default:
        throw new Error(UNSUPPORTED_TYPE_ERROR_MSG);
      }
    }
  }

  return table;
}

export function getInputOutputColumns(inputs: DG.Column[], outputs: DG.Column[]): DG.Column[] {
  outputs.forEach((outCol) => {
    inputs.forEach((inCol) => {
      if (inCol.name === outCol.name) {        
        inCol.name = `${inCol.name} (input)`;
        outCol.name = `${outCol.name} (output)`;
      }
    })
  });

  return inputs.concat(outputs);
};
