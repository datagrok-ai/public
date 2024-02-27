//

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const UNSUPPORTED_TYPE_ERROR_MSG = 'Unsupported output type.';
const FUNCCALL_ARRAY_IS_EMPTY_MSG = 'Funccall array is empty.';
const INCORRECT_ROW_MSG = 'Incorrect row number.';

type CellCoordinates = {
  idx: number,
  columnName: string,
};

export type OutputDataFromUI = {
  prop: DG.Property,
  value: {
    row: number | null,
    colName: string,
    colValue: number,
  },
};

export type OutputInfo = {
  prop: DG.Property,
  elements: CellCoordinates [],
};

export type SensitivityAnalysisResult = {
  funcEvalResults: DG.DataFrame,
  funcCalls: DG.FuncCall[],
};

function getCellCoordinatesFromRow(df: DG.DataFrame, spec: OutputDataFromUI): CellCoordinates[] {
  const rowCount = df.rowCount;
  const row = spec.value.row;

  // Last/first row case
  if (row !== null)
    return df.columns.names().map((colName) => ({idx: row, columnName: colName}));

  // By value in column
  const col = df.col(spec.value.colName);

  if (col === null) {
    // eslint-disable-next-line max-len
    grok.shell.warning(`Output dataframe "${df.name}" doesn't contain the "${spec.value.colName}" column. The last row is analyzed.`);
    return df.columns.names().map((colName) => ({idx: rowCount - 1, columnName: colName}));
  }

  if ((col.type !== DG.COLUMN_TYPE.INT) && (col.type !== DG.COLUMN_TYPE.FLOAT)) {
    // eslint-disable-next-line max-len
    grok.shell.warning(`Non-supported type of the "${spec.value.colName}" column: ${col.type}. Last row of the output dataframe "${df.name}" is analyzed.`);
    return df.columns.names().map((colName) => ({idx: rowCount - 1, columnName: colName}));
  }

  const raw = col.getRawData();
  const val = spec.value.colValue;
  let closestIdx = 0;
  let closestDist = Math.abs(raw[0] - val);
  let curDist: number;

  for (let i = 1; i < rowCount; ++i) {
    curDist = Math.abs(raw[i] - val);

    if (curDist < closestDist) {
      closestDist = curDist;
      closestIdx = i;
    }
  }

  return df.columns.names().map((colName) => ({idx: closestIdx, columnName: colName}));
} // getCellCoordinatesFromRow

function getExtendedOutputsSpecification(funcCall: DG.FuncCall,
  outputsSpecification: OutputDataFromUI[]): OutputInfo[] {
  const extendedOutputsSpecification = [] as OutputInfo[];

  for (const spec of outputsSpecification) {
    if (spec.prop.propertyType == DG.TYPE.DATA_FRAME) {
      extendedOutputsSpecification.push({
        prop: spec.prop,
        elements: getCellCoordinatesFromRow(funcCall.outputs[spec.prop.name], spec),
      });
    } else {
      extendedOutputsSpecification.push({
        prop: spec.prop,
        elements: [],
      });
    }
  }

  return extendedOutputsSpecification;
}

function getEmptyOutputTable(funcCall: DG.FuncCall,
  outputsSpecification: OutputInfo[], rowCount: number): DG.DataFrame {
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

export function getOutput(funcCalls: DG.FuncCall[], outputsSpecification: OutputDataFromUI[]): DG.DataFrame {
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
        const lastIdx = df.rowCount - 1;

        if (item.elements) {
          for (const elem of item.elements) {
            table.col(elem.columnName)?.set(
              row,
              df.col(elem.columnName)?.get((elem.idx !== -1) ? elem.idx : lastIdx),
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
    });
  });

  return inputs.concat(outputs);
};

export function getDataFrameFromInputsOutputs(inputs: DG.Column[], outputs: DG.Column[]): DG.DataFrame {
  const cols = getInputOutputColumns(inputs, outputs);
  const len = cols.length;

  if (len === 0)
    return DG.DataFrame.create();

  const df = DG.DataFrame.fromColumns([cols[0]]);

  for (let i = 1; i < len; ++i) {
    cols[i].name = df.columns.getUnusedName(cols[i].name);
    df.columns.add(cols[i]);
  }

  return df;
}
