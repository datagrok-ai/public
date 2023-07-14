import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {UiUtils} from '@datagrok-libraries/compute-utils';
import {CampaignFieldTypes, ITemplate, IngestType} from '../types';
import * as C from '../consts';
import '../../../css/hit-triage.css';

type INewCampaignResult = {
    df: DG.DataFrame,
    type: IngestType,
    campaignProps: {[key: string]: any}
}

type ICampaignAccordeon = {
    promise: Promise<INewCampaignResult>,
    root: HTMLElement,
    cancelPromise: Promise<void>
}
export function newCampaignAccordeon(template: ITemplate): ICampaignAccordeon {
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

  const fileInput = UiUtils.fileInput('', null, (f: File) => onFileChange(f), undefined);
  const fileInputDiv = ui.divV([fileInput, errorDiv]);
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

  const campaignProps = template.campaignFields.map((field) =>
    DG.Property.fromOptions({name: field.name, type: CampaignFieldTypes[field.type], nullable: !field.required}));

  const campaignPropsObject: {[key: string]: any} = {};
  const campaignPropsForm = ui.input.form(campaignPropsObject, campaignProps);
  if (template.dataSourceType === 'File') {
    fileInputDiv.style.display = 'block';
    functionInputDiv.style.display = 'none';
  } else {
    fileInputDiv.style.display = 'none';
    functionInputDiv.style.display = 'block';
  }

  const accordeon = ui.accordion();
  accordeon.root.classList.add('hit-triage-new-campaign-accordeon');
  accordeon.addPane('File source', () => dataInputsDiv, true);
  campaignProps.length && accordeon.addPane('Campaign info', () => campaignPropsForm, true);
  const content = ui.div(accordeon.root);
  const buttonsDiv = ui.divH([]);
  content.appendChild(buttonsDiv);
  const promise = new Promise<INewCampaignResult>((resolve) => {
    const onOkProxy = async () => {
      if (template.dataSourceType === 'File') {
        if (!file) {
          grok.shell.error('No file selected');
          return;
        }
        const df = DG.DataFrame.fromCsv(await file.text());
        df.name = file.name;
        resolve({df, type: 'File', campaignProps: campaignPropsObject});
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
        resolve({df, type: 'Query', campaignProps: campaignPropsObject});
      };
    };
    const startCampaignButton = ui.bigButton('Start campaign', () => onOkProxy());
    buttonsDiv.appendChild(startCampaignButton);
  });

  const cancelPromise = new Promise<void>((resolve) => {
    const cancelButton = ui.bigButton('Cancel', () => resolve());
    cancelButton.classList.add('hit-triage-accordeon-cancel-button');
    buttonsDiv.appendChild(cancelButton);
  });

  return {promise, root: content, cancelPromise};
}
