import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {CampaignFieldTypes, HitTriageCampaignField, HitTriageCampaignFieldType,
  IComputeDialogResult, HitTriageTemplate, IngestType, INewTemplateResult} from '../types';
import * as C from '../consts';
import '../../../css/hit-triage.css';
import {chemFunctionsDialog} from '../dialogs/functions-dialog';


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

  const ingestTypeInput = ui.choiceInput<IngestType>('Ingest using', 'Query', ['Query', 'File']);

  const fieldsEditor = getCampaignFieldEditors();

  const form = ui.divV([
    ui.h2('Details'),
    ui.divV([templateNameInput, errorDiv]),
    ui.divV([templateKeyInput, keyErrorDiv]),
    ingestTypeInput.root,
    fieldsEditor.fieldsDiv,
    ui.h2('Compute'),
    funcInput.root,
    ui.h2('Submit'),
    submitFunctionInput.root,
  ], 'ui-form');

  const content = ui.div(form);
  const buttonsDiv = ui.divH([]);
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
      };
      saveTemplate(out);
      grok.shell.info('Template created successfully');
      resolve(out);
    }
    const createTemplateButton = ui.bigButton(C.i18n.createTemplate, () => onOkProxy());
    buttonsDiv.appendChild(createTemplateButton);
  });

  const cancelPromise = new Promise<void>((resolve) => {
    const cancelButton = ui.button(C.i18n.cancel, () => resolve());
    cancelButton.classList.add('hit-triage-accordeon-cancel-button');
    //buttonsDiv.appendChild(cancelButton);
  });
  return {root: content, template: promise, cancelPromise};
}

export function getCampaignFieldEditors() {
  const getNewFieldEditor = () => {
    const nameInput = ui.stringInput('Name', '', () => out.changed = true);
    const typeInput = ui.choiceInput('Type', 'String', Object.keys(CampaignFieldTypes), () => out.changed = true);
    const requiredInput = ui.boolInput('Required', false, () => out.changed = true);
    requiredInput.classList.add('mx-5');
    const addFieldButton = ui.icons.add(() => {
      if (!fields.length || fields[fields.length - 1].changed) {
        const newField = getNewFieldEditor();
        fields.push(newField);
        fieldsContainer.appendChild(getFieldDiv(newField));
      }
    }, 'Add field');
    addFieldButton.classList.add('hit-triage-add-campaign-field-button');
    const out = {changed: false, nameInput, typeInput, requiredInput, addFieldButton};
    return out;
  };

  const fields: ReturnType<typeof getNewFieldEditor>[] = [getNewFieldEditor()];

  function getFieldParams(): HitTriageCampaignField[] {
    return fields.filter((f) => f.changed && f.nameInput.value).map((f) => ({
      name: f.nameInput.value,
      type: f.typeInput.value as HitTriageCampaignFieldType,
      required: f.requiredInput.value ?? false,
    }));
  };

  function getFieldDiv(field: ReturnType<typeof getNewFieldEditor>) {
    const removeButton = ui.icons.delete(() => {
      if (fields.length <= 1)
        return;
      fieldDiv.remove();
      fields.splice(fields.indexOf(field), 1);
    }, 'Remove Field');
    removeButton.classList.add('hit-triage-remove-campaign-field-button');
    field.requiredInput.addOptions(removeButton);
    field.requiredInput.addOptions(field.addFieldButton);
    const fieldDiv = ui.divH([
      field.nameInput.root,
      field.typeInput.root,
      field.requiredInput.root,
      // removeButton,
      // field.addFieldButton,
    ], {classes: 'hit-triage-campaign-field-div ui-input-row ui-input-root'});
    return fieldDiv;
  }
  const fieldsContainer = ui.divV([ui.divText('Additional Fields', {classes: 'hit-triage-additional-fields-title'}),
    getFieldDiv(fields[0])]);
  return {
    getFields: getFieldParams,
    fieldsDiv: ui.divV([fieldsContainer]),
  };
}

export function saveTemplate(template: HitTriageTemplate) {
  _package.files.writeAsText(`Hit Triage/templates/${template.name}.json`, JSON.stringify(template));
}
