import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import {saveAs} from 'file-saver';

export function exportFuncCall(call: DG.FuncCall) {
  //todo: check status
  // if (call.status != FuncCall.STATUS_COMPLETED) ...

  const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
  const exportWorkbook = new ExcelJS.Workbook();

  const isScalarType = (type: DG.TYPE) => {
    return DG.TYPES_SCALAR.has(type);
  };

  const isDataframe = (type: DG.TYPE) => (type === DG.TYPE.DATA_FRAME);

  const dfOutputs = call.func.outputs.filter(
    (output: DG.Property) => isDataframe(output.propertyType),
  ) as DG.Property[];

  const scalarOutputs = call.func.outputs.filter(
    (output: DG.Property) => isScalarType(output.propertyType),
  ) as DG.Property[];

  dfOutputs.forEach((dfOutput) => {
    const currentDfSheet = exportWorkbook.addWorksheet(dfOutput.name);
    const currentDf = (call.outputs[dfOutput.name] as DG.DataFrame);
    currentDfSheet.addRow((currentDf.columns as DG.ColumnList).names());
    for (let i = 0; i < currentDf.rowCount; i++) {
      currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: DG.Cell) => cell.value));
    }
  });

  const scalarsSheet = exportWorkbook.addWorksheet('Scalars');
  scalarOutputs.forEach((scalarOutput) => {
    scalarsSheet.addRow([scalarOutput.name, call.outputs[scalarOutput.name]]);
  });

  exportWorkbook.xlsx.writeBuffer().then((data) => {
    const blob = new Blob([data], {type: BLOB_TYPE});
    saveAs(blob, call.func.name);
  });
}
