/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model_handler';
import {selectOutliersManually} from './outliers_selection';
import {Cell, ColumnList, DataFrame, FuncCall, Property} from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import {saveAs} from 'file-saver';

export const _package = new DG.Package();

//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export function init() {
  console.log('init');
  DG.ObjectHandler.register(new ModelHandler());
}

//name: Model Catalog
//tags: app
export function modelCatalog() {
  const v = DG.CardView.create({dataSource: grok.dapi.scripts, permanentFilter: '#model'});
  v.meta = new ModelHandler();
  v.name = 'Models';
  v.permanentFilter = '#model';
  grok.shell.addView(v);
}

//name: manualOutlierDetectionDialog
//input: dataframe inputData
//output: dataframe augmentedInput
//output: dataframe editedInput
export async function manualOutlierSelectionDialog(inputData: DataFrame) {
  const {augmentedInput, editedInput} = await selectOutliersManually(inputData);
  return {augmentedInput, editedInput};
}


//name: testExportFuncCall
export async function testExportFuncCall() {
  const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
  const funcName = 'Compute:Chemometric';
  const workbook = new ExcelJS.Workbook();

  const isScalarType = (type: string) => {
    const SCALAR_TYPES = ['int', 'double', 'bool', 'string', 'datetime'];
    return SCALAR_TYPES.includes(type);
  };

  const isDataframe = (type: string) => (type === 'dataframe');

  grok.functions.eval(funcName).then((res) => {
    const targetFuncScript = res.prepare({a: 3, b: 3});
    const dfOutputs = res.outputs.filter((output: Property) => isDataframe(output.propertyType)) as Property[];
    const scalarOutputs = res.outputs.filter((output: Property) => isScalarType(output.propertyType)) as Property[];
    targetFuncScript.call().then((data: FuncCall) => {
      dfOutputs.forEach((dfOutput) => {
        const currentDfSheet = workbook.addWorksheet(dfOutput.name);
        const currentDf = (data.outputs[dfOutput.name] as DataFrame);
        currentDfSheet.addRow((currentDf.columns as ColumnList).names());
        for (let i =0; i< currentDf.rowCount; i++) {
          currentDfSheet.addRow([...currentDf.row(i).cells].map((cell: Cell) => cell.value));
        }
      });
      const scalarsSheet = workbook.addWorksheet('Scalars');
      scalarOutputs.forEach((scalarOutput) => {
        scalarsSheet.addRow([scalarOutput.name, data.outputs[scalarOutput.name]]);
      });
      workbook.xlsx.writeBuffer().then((data) => {
        const blob = new Blob([data], {type: BLOB_TYPE});
        saveAs(blob, funcName);
      });
    });
  });
}
