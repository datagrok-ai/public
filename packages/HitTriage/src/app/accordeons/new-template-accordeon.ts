import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {CampaignFieldTypes, ICampaignField, ICampaignFieldType, ITemplate, IngestType} from '../types';
import * as C from '../consts';
import '../../../css/hit-triage.css';

type INewTemplateResult = {
    template: Promise<ITemplate>,
    root: HTMLElement,
    cancelPromise: Promise<void>,
}

export async function createTemplateAccordeon(): Promise<INewTemplateResult> {
  const functions = DG.Func.find({tags: [C.HitTriageComputeFunctionTag]});
  const availableTemplates = (await _package.files.list('templates')).map((file) => file.name.slice(0, -5));
  const availableTemplateKeys: string[] = [];
  for (const tn of availableTemplates) {
    const t: ITemplate = JSON.parse(await _package.files.readAsText(`templates/${tn}.json`));
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

  templateKeyInput.setTooltip('Template key used for campeign prefix');
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
  const functionsInput = ui.multiChoiceInput('Select functions', [],
    Object.keys(functionsMap), () => {});
  functionsInput.setTooltip('Select functions to be applied to the data');
  functionsInput.root.classList.add('hit-triage-new-template-functions-input');

  const ingestTypeInput = ui.choiceInput<IngestType>('Ingest using', 'Query', ['Query', 'File']);

  const fieldsEditor = getCampaignFieldEditors();

  const detailsDiv = ui.divV(
    [ui.divV([templateNameInput, errorDiv]), ui.divV([templateKeyInput, keyErrorDiv]), fieldsEditor.fieldsDiv]);
  const accordeon = ui.accordion();
  accordeon.root.classList.add('hit-triage-new-template-accordeon');
  accordeon.addPane('Details', () => detailsDiv, true);
  accordeon.addPane('Ingest', () => ingestTypeInput.root, true);
  accordeon.addPane('Compute', () => functionsInput.root, true);
  accordeon.addPane('Submit', () => submitFunctionInput.root, true);

  const content = ui.div(accordeon.root);
  const buttonsDiv = ui.divH([]);
  content.appendChild(buttonsDiv);
  const promise = new Promise<ITemplate>((resolve) => {
    async function onOkProxy() {
      if (errorDiv.style.opacity === '100%') {
        grok.shell.error('Template name is empty or already exists');
        return;
      }
      const submitFunction = submitFunctionInput.value ? submitFunctionsMap[submitFunctionInput.value] : undefined;
      const out: ITemplate = {
        name: templateNameInput.value,
        key: templateKeyInput.value,
        campaignFields: fieldsEditor.getFields(),
        dataSourceType: ingestTypeInput.value ?? 'Query',
        compute: {
          descriptors: {
            enabled: functionsInput.value!.includes('Descriptors'),
          },
          functions: functionsInput.value!.filter((f) => f !== 'Descriptors')
            .map((f) => {
              const funcInfo = functionsMap[f].split(':');
              return {
                package: funcInfo[0],
                name: funcInfo[1],
              };
            }),
        },
        ...(submitFunction ? {submit: {fName: submitFunction.name, package: submitFunction.package.name}} : {}),
      };
      saveTemplate(out);
      resolve(out);
    }
    const createTemplateButton = ui.bigButton('Create template', () => onOkProxy());
    buttonsDiv.appendChild(createTemplateButton);
  });

  const cancelPromise = new Promise<void>((resolve) => {
    const cancelButton = ui.bigButton('Cancel', () => resolve());
    cancelButton.classList.add('hit-triage-accordeon-cancel-button');
    buttonsDiv.appendChild(cancelButton);
  });
  return {root: content, template: promise, cancelPromise};
}

export function getCampaignFieldEditors() {
  const getNewFieldEditor = () => {
    const nameInput = ui.stringInput('Name', '', () => out.changed = true);
    const typeInput = ui.choiceInput('Type', 'String', Object.keys(CampaignFieldTypes), () => out.changed = true);
    const requiredInput = ui.boolInput('Required', false, () => out.changed = true);
    requiredInput.classList.add('mx-5');
    const out = {changed: false, nameInput, typeInput, requiredInput};
    return out;
  };

  const fields: ReturnType<typeof getNewFieldEditor>[] = [getNewFieldEditor()];

  function getFieldParams(): ICampaignField[] {
    return fields.filter((f) => f.changed && f.nameInput.value).map((f) => ({
      name: f.nameInput.value,
      type: f.typeInput.value as ICampaignFieldType,
      required: f.requiredInput.value ?? false,
    }));
  };

  function getFieldDiv(field: ReturnType<typeof getNewFieldEditor>) {
    const fieldDiv = ui.divH([
      field.nameInput.root,
      field.typeInput.root,
      field.requiredInput.root,
      ui.button('X', () => {
        fieldDiv.remove();
        fields.splice(fields.indexOf(field), 1);
      }, 'Remove Field'),
    ], {classes: 'hit-triage-campaign-field-div'});
    return fieldDiv;
  }
  const fieldsContainer = ui.divV([ui.divText('Additional Fields', {classes: 'hit-triage-additional-fields-title'}),
    getFieldDiv(fields[0])]);
  const addFieldButton = ui.bigButton('+', () => {
    if (!fields.length || fields[fields.length - 1].changed) {
      const newField = getNewFieldEditor();
      fields.push(newField);
      fieldsContainer.appendChild(getFieldDiv(newField));
    }
  }, 'Add field');
  addFieldButton.classList.add('hit-triage-add-campaign-field-button');
  return {
    getFields: getFieldParams,
    fieldsDiv: ui.divV([fieldsContainer, addFieldButton]),
  };
}


export function saveTemplate(template: ITemplate) {
  _package.files.writeAsText(`templates/${template.name}.json`, JSON.stringify(template));
}
