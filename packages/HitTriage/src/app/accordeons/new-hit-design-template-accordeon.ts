import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import * as C from '../consts';
import '../../../css/hit-triage.css';
import {HitDesignTemplate, IComputeDialogResult, INewTemplateResult} from '../types';
import {chemFunctionsDialog} from '../dialogs/functions-dialog';
import {getCampaignFieldEditors} from './new-template-accordeon';
import {ItemType, ItemsGrid} from '@datagrok-libraries/utils/src/items-grid';

export async function newHitDesignTemplateAccordeon(): Promise<INewTemplateResult<HitDesignTemplate>> {
  const functions = DG.Func.find({tags: [C.HitTriageComputeFunctionTag]});
  const availableTemplates = (await _package.files.list('Hit Design/templates'))
    .map((file) => file.name.slice(0, -5)); // Remove .json from end

  const availableTemplateKeys: string[] = [];
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
  // ######### TEMPLATE NAME INPUT #########
  const templateNameInput = ui.stringInput('Name', '', () => {
    if (templateNameInput.value === '' || availableTemplates.includes(templateNameInput.value)) {
      templateNameInput.root.style.borderBottom = '1px solid red';
      errorDiv.style.opacity = '100%';
    } else {
      templateNameInput.root.style.borderBottom = 'none';
      errorDiv.style.opacity = '0%';
    }
  });
  // ######### TEMPLATE KEY INPUT #########
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
  } as unknown as HitDesignTemplate;
  const funcInput = await chemFunctionsDialog((res) => {funcDialogRes = res;}, () => null,
    dummyTemplate, false);
  funcInput.root.classList.add('hit-triage-new-template-functions-input');
  const fieldsEditor = getCampaignFieldEditors();
  const tileCategoriesEditor = getTileCategoryEditor();
  const detailsDiv = ui.divV(
    [ui.divV([templateNameInput, errorDiv]), ui.divV([templateKeyInput, keyErrorDiv]), fieldsEditor.fieldsDiv]);

  const form = ui.div(
    [ui.h2('Details'),
      detailsDiv,
      ui.h2('Stages'),
      tileCategoriesEditor.fieldsDiv,
      ui.h2('Compute'),
      funcInput.root,
      ui.h2('Submit'),
      submitFunctionInput.root,
    ], 'ui-form');
  const buttonsDiv = ui.buttonsInput([]);
  form.appendChild(buttonsDiv);
  const buttonsContainerDiv = buttonsDiv.getElementsByClassName('ui-input-editor')?.[0] ?? buttonsDiv;
  const cancelPromise = new Promise<void>((resolve) => {
    const cancelButton = ui.button(C.i18n.cancel, () => resolve());
    buttonsContainerDiv.appendChild(cancelButton);
  });
  const promise = new Promise<HitDesignTemplate>((resolve) => {
    async function onOkProxy() {
      funcInput.okProxy();
      if (errorDiv.style.opacity === '100%' || !templateNameInput.value || templateNameInput.value === '') {
        grok.shell.error('Template name is empty or already exists');
        return;
      }

      if (keyErrorDiv.style.opacity === '100%' || !templateKeyInput.value || templateKeyInput.value === '') {
        grok.shell.error('Template key is empty or already exists');
        return;
      }

      const submitFunction = submitFunctionInput.value ? submitFunctionsMap[submitFunctionInput.value] : undefined;
      const out: HitDesignTemplate = {
        name: templateNameInput.value,
        key: templateKeyInput.value,
        campaignFields: fieldsEditor.getFields(),
        stages: tileCategoriesEditor.getFields(),
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
      saveHitDesignTemplate(out);
      grok.shell.info('Template created successfully');
      resolve(out);
    }
    const createTemplateButton = ui.bigButton(C.i18n.createTemplate, () => onOkProxy());
    buttonsContainerDiv.appendChild(createTemplateButton);
  });

  return {root: form, template: promise, cancelPromise};
}

export function saveHitDesignTemplate(template: HitDesignTemplate) {
  _package.files.writeAsText(`Hit Design/templates/${template.name}.json`, JSON.stringify(template));
}


export function getTileCategoryEditor() {
  const props = [DG.Property.fromOptions({name: 'Name', type: DG.TYPE.STRING})];
  const itemsGrid = new ItemsGrid(props, undefined, {horizontalInputNames: true});
  let addingItem: ItemType = {};
  function getFieldParams(): string[] {
    const items = itemsGrid.items.filter((f) => f.Name).map((f) => f.Name);
    if (addingItem.Name && addingItem.Name !== '')
      items.push(addingItem.Name);
    return items;
  }
  itemsGrid.onItemAdded.subscribe((item) => {
    addingItem = item ?? {};
  });
  itemsGrid.onAddingItemChanged.subscribe((item) => {
    if (item)
      addingItem = item.item;
  });

  return {
    getFields: getFieldParams,
    fieldsDiv: itemsGrid.root,
  };
}
