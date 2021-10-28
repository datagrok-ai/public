import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Cell, ColumnList, DataFrame, FuncCall, Property, TYPE, TYPES_SCALAR} from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import {saveAs} from 'file-saver';

export function exportFuncCall(targetFuncScript: FuncCall) {
  const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
  const exportWorkbook = new ExcelJS.Workbook();

  const isScalarType = (type: TYPE) => {
    return TYPES_SCALAR.has(type);
  };

  const isDataframe = (type: TYPE) => (type === TYPE.DATA_FRAME);

  const dfOutputs = targetFuncScript.func.outputs.filter(
    (output: Property) => isDataframe(output.propertyType),
  ) as Property[];
  const scalarOutputs = targetFuncScript.func.outputs.filter(
    (output: Property) => isScalarType(output.propertyType),
  ) as Property[];

  targetFuncScript.call().then((data: FuncCall) => {
    dfOutputs.forEach((dfOutput) => {
      const currentDfSheet = exportWorkbook.addWorksheet(dfOutput.name);
      const currentDf = (data.outputs[dfOutput.name] as DataFrame);
      currentDfSheet.addRow((currentDf.columns as ColumnList).names());
      for (let i =0; i< currentDf.rowCount; i++) {
        currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: Cell) => cell.value));
      }
    });

    const scalarsSheet = exportWorkbook.addWorksheet('Scalars');
    scalarOutputs.forEach((scalarOutput) => {
      scalarsSheet.addRow([scalarOutput.name, data.outputs[scalarOutput.name]]);
    });

    exportWorkbook.xlsx.writeBuffer().then((data) => {
      const blob = new Blob([data], {type: BLOB_TYPE});
      saveAs(blob, targetFuncScript.func.name);
    });
  });
}
