import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {CampaignFieldTypes, HitTriageCampaignField, HitTriageCampaignFieldType,
  IComputeDialogResult, HitTriageTemplate, IngestType, INewTemplateResult} from '../types';
import * as C from '../consts';
import '../../../css/hit-triage.css';
import {chemFunctionsDialog} from '../dialogs/functions-dialog';
import {ItemsGrid} from '@datagrok-libraries/utils/src/items-grid';


export async function createTemplateAccordeon(): Promise<INewTemplateResult<HitTriageTemplate>> {
  const functions = DG.Func.find({tags: [C.HitTriageComputeFunctionTag]});
  const availableTemplates = (await _package.files.list('Hit Triage/templates')).map((file) => file.name.slice(0, -5));
  const availableTemplateKeys: string[] = [];
  for (const tn of availableTemplates) {
    const t: HitTriageTemplate = JSON.parse(await _package.files.readAsText(`Hit Triage/templates/${tn}.json`));
    availableTemplateKeys.push(t.key);
  }

  const availableSubmitFunctions = DG.Func.find({tags: [C.HitTriageSubmitTag]});
  const submitFunctionsMap: {[key: string]: DG.Func} = {};
  availableSubmitFunctions.forEach((func) => {
    submitFunctionsMap[func.friendlyName ?? func.name] = func;
  });
  const submitFunctionInput = ui.choiceInput('Submit function', null, Object.keys(submitFunctionsMap));
  submitFunctionInput.value = null;
  submitFunctionInput.setTooltip('Select function to be called upon submitting');
  const errorDiv = ui.divText('Template name is empty or already exists', {classes: 'hit-triage-error-div'});

  const keyErrorDiv = ui.divText('Template key is empty or already exists', {classes: 'hit-triage-error-div'});

  const templateNameInput = ui.stringInput('Name', '', () => {
    if (templateNameInput.value === '' || availableTemplates.includes(templateNameInput.value)) {
      templateNameInput.root.style.borderBottom = '1px solid red';
      errorDiv.style.opacity = '100%';
    } else {
      templateNameInput.root.style.borderBottom = 'none';
      errorDiv.style.opacity = '0%';
    }
  });
  const templateKeyInput = ui.stringInput('Key', '', () => {
    if (templateKeyInput.value === '' || availableTemplateKeys.includes(templateKeyInput.value)) {
      templateKeyInput.root.style.borderBottom = '1px solid red';
      keyErrorDiv.style.opacity = '100%';
    } else {
      templateKeyInput.root.style.borderBottom = 'none';
      keyErrorDiv.style.opacity = '0%';
    }
  });

  templateKeyInput.setTooltip('Template key used for campaign prefix');
  templateNameInput.setTooltip('Template name');
  templateNameInput.root.style.borderBottom = 'none';
  templateKeyInput.root.style.borderBottom = 'none';
  errorDiv.style.opacity = '0%';
  keyErrorDiv.style.opacity = '0%';

  const functionsMap: {[key: string]: string} = {};
  functionsMap['Descriptors'] = 'Descriptors';
  functions.forEach((func) => {
    functionsMap[func.friendlyName ?? func.name] = `${func.package.name}:${func.name}`;
  });

  let funcDialogRes: IComputeDialogResult | null = null;
  // used just for functions editor
  const dummyTemplate = {
    compute: {
      descriptors: {
        enabled: true,
        args: [],
      },
      functions: functions.map((f) => ({name: f.name, package: f.package.name, args: []})),
    },
  } as unknown as HitTriageTemplate;
  const funcInput = await chemFunctionsDialog((res) => {funcDialogRes = res;}, () => null,
    dummyTemplate, false);
  funcInput.root.classList.add('hit-triage-new-template-functions-input');
  // functions that have special tag and are applicable for data source. they should return a dataframe with molecules
  const dataSourceFunctions = DG.Func.find({tags: [C.HitTriageDataSourceTag]});
  const dataSourceFunctionsMap: {[key: string]: DG.Func} = {};
  dataSourceFunctions.forEach((func) => {
    dataSourceFunctionsMap[func.friendlyName ?? func.name] = func;
  });
  const dataSourceFunctionInput = ui.choiceInput(
    C.i18n.dataSourceFunction, Object.keys(dataSourceFunctionsMap)[0],
    Object.keys(dataSourceFunctionsMap));
  const ingestTypeInput = ui.choiceInput<IngestType>('Ingest using', 'Query', ['Query', 'File'], () => {
    if (ingestTypeInput.value !== 'Query')
      dataSourceFunctionInput.root.style.display = 'none';
    else
      dataSourceFunctionInput.root.style.display = 'block';
  });

  const fieldsEditor = getCampaignFieldEditors();

  const form = ui.divV([
    ui.h3('Details'),
    ui.divV([templateNameInput, errorDiv]),
    ui.divV([templateKeyInput, keyErrorDiv]),
    ingestTypeInput.root,
    dataSourceFunctionInput.root,
    fieldsEditor.fieldsDiv,
    ui.h3('Compute'),
    funcInput.root,
    ui.h3('Submit'),
    submitFunctionInput.root,
  ], 'ui-form');

  const buttonsDiv = ui.buttonsInput([]);
  const buttonsContainerDiv = buttonsDiv.getElementsByClassName('ui-input-editor')?.[0] ?? buttonsDiv;
  form.appendChild(buttonsDiv);
  const promise = new Promise<HitTriageTemplate>((resolve) => {
    async function onOkProxy() {
      funcInput.okProxy();
      if (errorDiv.style.opacity === '100%') {
        grok.shell.error('Template name is empty or already exists');
        return;
      }
      const submitFunction = submitFunctionInput.value ? submitFunctionsMap[submitFunctionInput.value] : undefined;
      const out: HitTriageTemplate = {
        name: templateNameInput.value,
        key: templateKeyInput.value,
        campaignFields: fieldsEditor.getFields(),
        dataSourceType: ingestTypeInput.value ?? 'Query',
        compute: {
          descriptors: {
            enabled: !!funcDialogRes?.descriptors?.length,
            args: funcDialogRes?.descriptors ?? [],
          },
          functions: Object.entries(funcDialogRes?.externals ?? {}).map(([funcName, args]) => {
            const splitFunc = funcName.split(':');
            return ({
              name: splitFunc[1],
              package: splitFunc[0],
              args: args,
            });
          }),
        },
        ...(submitFunction ? {submit: {fName: submitFunction.name, package: submitFunction.package.name}} : {}),
        queryFunctionName: (ingestTypeInput.value === 'Query') ? dataSourceFunctionInput.value ?? undefined : undefined,
      };
      saveTemplate(out);
      grok.shell.info('Template created successfully');
      resolve(out);
    }
    const createTemplateButton = ui.bigButton(C.i18n.createTemplate, () => onOkProxy());
    buttonsContainerDiv.appendChild(createTemplateButton);
  });

  const cancelPromise = new Promise<void>((resolve) => {
    const cancelButton = ui.button(C.i18n.cancel, () => resolve());
    buttonsContainerDiv.appendChild(cancelButton);
  });
  return {root: form, template: promise, cancelPromise};
}

export function getCampaignFieldEditors() {
  const props = [DG.Property.fromOptions({name: 'Name', type: DG.TYPE.STRING}),
    DG.Property.fromOptions({name: 'Type', type: DG.TYPE.STRING, choices: Object.keys(CampaignFieldTypes)}),
    DG.Property.fromOptions({name: 'Required', type: DG.TYPE.BOOL})];
  const itemsGrid = new ItemsGrid(props, undefined, {horizontalInputNames: true});

  function getFieldParams(): HitTriageCampaignField[] {
    return itemsGrid.items.filter((f) => f.Name).map((f) => ({
      name: f.Name,
      type: f.Type as HitTriageCampaignFieldType,
      required: f.Required ?? false,
    }));
  }
  const fieldsContainer = ui.divV([ui.h3('Additional fields'),
    itemsGrid.root]);
  return {
    getFields: getFieldParams,
    fieldsDiv: ui.divV([fieldsContainer]),
  };
}

export function saveTemplate(template: HitTriageTemplate) {
  _package.files.writeAsText(`Hit Triage/templates/${template.name}.json`, JSON.stringify(template));
}
