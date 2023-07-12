import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {UiUtils} from '@datagrok-libraries/compute-utils';
import {IngestType} from '../types';
import * as C from '../consts';

type IFileInputDialogResult = {
    df: DG.DataFrame,
    type: IngestType,
}
export function fileInputDialog(): Promise<IFileInputDialogResult> {
  return new Promise<IFileInputDialogResult>(async (resolve) => {
    let file: File | null = null;
    const errorDiv = ui.divText('', {style: {color: 'red'}});
    const onFileChange = async (f: File) => {
      try {
        const df = DG.DataFrame.fromCsv(await f.text());
        await df.meta.detectSemanticTypes();
        const molcol = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
        if (!molcol) {
          errorDiv.innerText = 'No molecules column found';
          return;
        }
        file = f;
        errorDiv.innerText = '';
      } catch (e) {
        errorDiv.innerText = 'Error parsing file';
      }
    };

    const onFileSourceChange = (source: IngestType) => {
      if (source === 'File') {
        fileInputDiv.style.display = 'block';
        functionInputDiv.style.display = 'none';
      } else {
        fileInputDiv.style.display = 'none';
        functionInputDiv.style.display = 'block';
      }
    };
    const fileInput = UiUtils.fileInput('', null, (f: File) => onFileChange(f), undefined);
    const fileInputDiv = ui.divV([fileInput, errorDiv]);
    const dataSourceInput = ui.choiceInput<IngestType>('Data source', 'File', ['File', 'Function'], onFileSourceChange);
    const dataSourceFunctions = DG.Func.find({tags: [C.HitTriageDataSourceTag]});
    const dataSourceFunctionsMap: {[key: string]: DG.Func} = {};
    dataSourceFunctions.forEach((func) => {
      dataSourceFunctionsMap[func.friendlyName ?? func.name] = func;
    });

    let funcCall: DG.FuncCall | null = null;
    const onDataFunctionChange = async () => {
      const func = dataSourceFunctionsMap[dataSourceFunctionInput.value!];
      funcCall = func.prepare();
      const editor = await funcCall.getEditor();
      funcEditorDiv.innerHTML = '';
      funcEditorDiv.appendChild(editor);
    };
    const funcEditorDiv = ui.div();
    const dataSourceFunctionInput = ui.choiceInput(
      'Data source function', Object.keys(dataSourceFunctionsMap)[0],
      Object.keys(dataSourceFunctionsMap), onDataFunctionChange);
    onDataFunctionChange();
    const functionInputDiv = ui.divV([dataSourceFunctionInput, funcEditorDiv]);
    functionInputDiv.style.display = 'none';
    const dataInputsDiv = ui.divV([fileInputDiv, functionInputDiv]);
    const content = ui.divV([dataSourceInput.root, dataInputsDiv]);

    const onOkProxy = async () => {
      if (dataSourceInput.value === 'File') {
        if (!file) {
          grok.shell.error('No file selected');
          return;
        }
        const df = DG.DataFrame.fromCsv(await file.text());
        df.name = file.name;
        resolve({df, type: 'File'});
      } else {
        const func = dataSourceFunctionsMap[dataSourceFunctionInput.value!];
        if (!func) {
          grok.shell.error('No function selected');
          return;
        }
        const funcCallInputs: {[key: string]: any} = {};
        Object.entries(funcCall!.inputs).forEach(([key, value]) => {
          funcCallInputs[key] = value;
        });
        const df: DG.DataFrame = await func.apply(funcCallInputs);
        resolve({df, type: 'Function'});
      };
    };

    ui.dialog('Select data source')
      .add(content)
      .onOK(onOkProxy)
      .show();
  });
}
