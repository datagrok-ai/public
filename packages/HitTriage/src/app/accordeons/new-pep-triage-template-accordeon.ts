/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {_package} from '../../package';
import {IComputeDialogResult, PepTriageTemplate, IngestType, INewTemplateResult} from '../types';
import * as C from '../consts';
import '../../../css/hit-triage.css';
import {chemFunctionsDialog} from '../dialogs/functions-dialog';
import {HitAppBase} from '../hit-app-base';
import {getLayoutInput} from './layout-input';
import {getFuncPackageNameSafe, loadTemplate} from '../utils';
import {getCampaignFieldEditors, showTemplateInputError} from './new-template-accordeon';
import {freeKey, freeName} from '../shared/template-editor';


export async function createPepTriageTemplateAccordeon(app: HitAppBase<any>,
  dataSourceFunctionMap: { [key: string]: DG.Func | DG.DataQuery | DG.Script },
  preset?: PepTriageTemplate,
): Promise<INewTemplateResult<PepTriageTemplate>> {
  const templatesFolder = `${app.appName}/templates`;
  const availableTemplates = (await _package.files.list(templatesFolder))
    .filter((file) => file.name.endsWith('.json'))
    .map((file) => file.name.slice(0, -5));
  const availableTemplateKeys: string[] = [];
  for (const tn of availableTemplates) {
    const t: PepTriageTemplate = await loadTemplate<PepTriageTemplate>(`${templatesFolder}/${tn}.json`);
    availableTemplateKeys.push(t.key);
  }

  const availableSubmitFunctions = Array.from(new Set(DG.Func.find({meta: {role: C.HitTriageSubmitTag}}).concat(DG.Func.find({tags: [C.HitTriageSubmitTag]}))));
  const submitFunctionsMap: {[key: string]: DG.Func} = {};
  availableSubmitFunctions.forEach((func) => {
    submitFunctionsMap[func.friendlyName ?? func.name] = func;
  });
  const presetSubmitKey = preset?.submit ? Object.keys(submitFunctionsMap)
    .find((k) => submitFunctionsMap[k]?.name === preset.submit?.fName) : null;
  const submitFunctionInput =
    ui.input.choice('Submit function', {value: presetSubmitKey ?? null, items: [null, ...Object.keys(submitFunctionsMap)]});
  if (!presetSubmitKey)
    submitFunctionInput.value = null;
  submitFunctionInput.nullable = true;
  submitFunctionInput.fireChanged();
  submitFunctionInput.setTooltip('Select function to be called upon submitting');
  const layoutInput = getLayoutInput(preset);
  const errorDiv = ui.divText('Template name is empty or already exists', {classes: 'hit-triage-error-div'});
  const keyErrorDiv = ui.divText('Template key is empty or already exists', {classes: 'hit-triage-error-div'});

  const nameProblem = () => !templateNameInput.value?.trim() || availableTemplates.includes(templateNameInput.value) ?
    'Template name is empty or already exists' : null;
  const keyProblem = () => !templateKeyInput.value?.trim() || availableTemplateKeys.includes(templateKeyInput.value) ?
    'Template key is empty or already exists' : null;
  const seqColProblem = () => seqColInput.value?.trim() ? null : 'Sequence column name is required';
  let createTemplateButton: HTMLButtonElement | null = null;
  let createProblem: string | null = null;
  let createDisabled = false;
  function updateCreate() {
    createProblem = nameProblem() ?? keyProblem() ?? seqColProblem();
    if (createTemplateButton && (createProblem !== null) !== createDisabled) {
      createDisabled = createProblem !== null;
      ui.setDisabled(createTemplateButton, createDisabled, () => createProblem);
    }
  }

  const templateNameInput = ui.input.string('Name', {value: preset ? freeName(preset.name, new Set(availableTemplates.map((n) => n.toLowerCase()))) : '', onValueChanged: () => {
    showTemplateInputError(templateNameInput, errorDiv, nameProblem() !== null);
    updateCreate();
  }});
  const templateKeyInput = ui.input.string('Key', {value: preset ? freeKey(preset.key, new Set(availableTemplateKeys.map((k) => k.toLowerCase()))) : '', onValueChanged: () => {
    showTemplateInputError(templateKeyInput, keyErrorDiv, keyProblem() !== null);
    updateCreate();
  }});

  templateKeyInput.setTooltip('Template key used for campaign prefix');
  templateNameInput.setTooltip('Template name');
  errorDiv.style.opacity = '0%';
  keyErrorDiv.style.opacity = '0%';

  // PepTriage-specific inputs
  const seqColInput = ui.input.string('Sequence Column', {value: preset?.sequenceColumnName ?? '', nullable: false});
  seqColInput.setTooltip('Name of the sequence column in the dataset (mandatory)');
  const seqColErrorDiv = ui.divText('Sequence column name is required', {classes: 'hit-triage-error-div'});
  seqColErrorDiv.style.opacity = '0%';
  seqColInput.onChanged.subscribe(() => {
    showTemplateInputError(seqColInput, seqColErrorDiv, seqColProblem() !== null);
    updateCreate();
  });

  const molColInput = ui.input.string('Molecule Column', {value: preset?.moleculeColumnName ?? '', nullable: true});
  molColInput.setTooltip('Name of the molecule column (optional). If empty, sequences will be auto-converted to molecules.');

  let funcDialogRes: IComputeDialogResult | null = null;
  const dummyTemplate = preset ? preset : {
    compute: {
      descriptors: {enabled: true, args: []},
      functions: [],
    },
  } as unknown as PepTriageTemplate;
  const funcInput = await chemFunctionsDialog(app, (res) => {funcDialogRes = res;}, () => null,
    dummyTemplate, false);
  funcInput.root.classList.add('hit-triage-new-template-functions-input');

  if (Object.entries(dataSourceFunctionMap).length === 0) {
    const dataSourceFunctions = Array.from(new Set(
      DG.Func.find({meta: {role: C.PepTriageDataSourceTag}}).concat(DG.Func.find({tags: [C.PepTriageDataSourceTag]}))));
    dataSourceFunctions.forEach((func) => {
      dataSourceFunctionMap[func.friendlyName ?? func.name] = func;
    });
  }
  const combinedSourceNames = Object.keys(dataSourceFunctionMap);
  const dataSourceFunctionInput = ui.input.choice(
    C.i18n.dataSourceFunction, {value: combinedSourceNames[0], items: combinedSourceNames});
  const ingestTypeInput = ui.input.choice<IngestType>('Ingest using', {value: preset?.dataSourceType ?? 'Query', items: ['Query', 'File'],
    onValueChanged: (value) => {
      dataSourceFunctionInput.root.style.display = value === 'Query' ? 'block' : 'none';
    }});
  if (preset?.dataSourceType)
    dataSourceFunctionInput.root.style.display = preset.dataSourceType === 'Query' ? 'block' : 'none';

  const fieldsEditor = getCampaignFieldEditors(preset?.campaignFields);

  const form = ui.div([
    ui.h2('Details'),
    ui.div([templateNameInput, errorDiv]),
    ui.div([templateKeyInput, keyErrorDiv]),
    ui.h2('Sequence Configuration'),
    ui.div([seqColInput, seqColErrorDiv]),
    molColInput.root,
    ui.h2('Data Ingestion'),
    ingestTypeInput.root,
    dataSourceFunctionInput.root,
    layoutInput.dataFileInput,
    fieldsEditor.fieldsDiv,
    ui.h2('Compute'),
    funcInput.root,
    ui.h2('Submit'),
    submitFunctionInput.root,
  ], 'ui-form');

  const buttonsDiv = ui.buttonsInput([]);
  const buttonsContainerDiv = buttonsDiv.getElementsByClassName('ui-input-editor')?.[0] ?? buttonsDiv;
  form.appendChild(buttonsDiv);
  const cancelPromise = new Promise<void>((resolve) => {
    const cancelButton = ui.button(C.i18n.cancel, () => resolve());
    buttonsContainerDiv.appendChild(cancelButton);
  });

  const promise = new Promise<PepTriageTemplate>((resolve) => {
    async function onOkProxy() {
      funcInput.okProxy();
      const problem = nameProblem() ?? keyProblem() ?? seqColProblem();
      if (problem) {
        grok.shell.error(problem);
        return;
      }
      const submitFunction = submitFunctionInput.value ? submitFunctionsMap[submitFunctionInput.value] : undefined;
      const out: PepTriageTemplate = {
        name: templateNameInput.value,
        key: templateKeyInput.value,
        sequenceColumnName: seqColInput.value!.trim(),
        moleculeColumnName: molColInput.value?.trim() || undefined,
        campaignFields: fieldsEditor.getFields(),
        dataSourceType: ingestTypeInput.value ?? 'Query',
        layoutViewState: layoutInput.getLayoutViewState() ?? undefined,
        compute: {
          descriptors: {
            enabled: !!funcDialogRes?.descriptors?.length,
            args: funcDialogRes?.descriptors ?? [],
          },
          functions: Object.entries(funcDialogRes?.externals ?? {}).map(([funcName, args]) => {
            const splitFunc = funcName.split(':');
            return ({name: splitFunc[1], package: splitFunc[0], args: args});
          }),
          scripts: Object.entries(funcDialogRes?.scripts ?? {})
            .filter(([name, _]) => name.startsWith(C.HTScriptPrefix) && name.split(':').length === 3)
            .map(([scriptId, args]) => {
              const scriptNameParts = scriptId.split(':');
              return ({name: scriptNameParts[1] ?? '', id: scriptNameParts[2] ?? '', args: args});
            }),
          queries: Object.entries(funcDialogRes?.queries ?? {})
            .filter(([name, _]) => name.startsWith(C.HTQueryPrefix) && name.split(':').length === 3)
            .map(([queryName, args]) => {
              const queryNameParts = queryName.split(':');
              return ({name: queryNameParts[1] ?? '', id: queryNameParts[2] ?? '', args: args});
            }),
        },
        ...(submitFunction ? {submit: {fName: submitFunction.name, package: getFuncPackageNameSafe(submitFunction)}} : {}),
        queryFunctionName: (ingestTypeInput.value === 'Query') ? dataSourceFunctionInput.value ?? undefined : undefined,
      };
      try {
        await _package.files.writeAsText(`${templatesFolder}/${out.name}.json`, JSON.stringify(out));
      }
      catch (e) {
        grok.shell.error(`Could not save template ${out.name}: ${e instanceof Error ? e.message : e}`);
        return;
      }
      grok.shell.info('Template created successfully');
      resolve(out);
    }
    createTemplateButton = ui.bigButton(C.i18n.createTemplate, () => onOkProxy());
    buttonsContainerDiv.appendChild(createTemplateButton);
    updateCreate();
  });

  return {root: form, template: promise, cancelPromise};
}
