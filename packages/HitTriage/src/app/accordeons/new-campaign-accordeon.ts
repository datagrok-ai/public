/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Scope, button, computed, signal} from '@datagrok-libraries/u2';
import type {LabelColumn} from '@datagrok-libraries/u2';
import {CampaignFieldTypes, HitTriageTemplate, IngestType} from '../types';
import * as C from '../consts';
import {HitAppBase} from '../hit-app-base';
import '../../../css/hit-triage.css';

export type INewCampaignResult = {
    df: DG.DataFrame,
    type: IngestType,
    campaignProps: {[key: string]: any},
    friendlyName?: string,
    extraData?: {[key: string]: any},
}

export type CampaignAccordeonOptions = {
    /** The scope the form lives in: it is built after an await, with no builder around. */
    owner: Scope,
    /** The label column of the Template row above the form: the form's captions set its width. */
    column?: LabelColumn,
    /** What Start does with the data; Start stays busy until it is done. */
    start: (result: INewCampaignResult) => Promise<unknown>,
    extraFormElements?: HTMLElement[],
    getExtraData?: () => {[key: string]: any},
}

const EMPTY = 'Value can\'t be empty';

/** Creates the new campaign form of Hit Triage and PepTriage: the campaign name, the data source (a file, or the
 * template's function and its parameters), the template's fields, and Start — disabled while the data or a required
 * value is missing, the reason its tooltip.
 * @param {HitTriageTemplate} template template for the campaign. it contains file source type and
 * additional campaign fields
 * @return {Promise<HTMLElement>} the form
 */
export async function newCampaignAccordeon(template: HitTriageTemplate,
  dataSourceFunctionsMap: {[key: string]: DG.Func | DG.DataQuery},
  dataSourceTag: string,
  options: CampaignAccordeonOptions): Promise<HTMLElement> {
  // the platform inputs report their edits here, and Start's reason is read again
  const edits = signal(0);
  const edited = () => edits.value++;
  const errorDiv = ui.divText('', {style: {color: 'red'}});
  // handling file input
  let fileDf: DG.DataFrame | null = null;
  const onFileChange = async () => {
    try {
      fileDf = dfInput.value;
      if (!fileDf)
        return;
      await fileDf.meta.detectSemanticTypes();
      const molcol = HitAppBase.moleculeColumn(fileDf);
      if (!molcol) {
        errorDiv.innerText = 'No molecules column found';
        return;
      }

      errorDiv.innerText = '';
    } catch (e) {
      errorDiv.innerText = 'Error parsing file';
    } finally {
      edited();
    }
  };

  const dfInput = ui.input.table('Dataframe', {onValueChanged: onFileChange});
  await onFileChange();
  const fileInputDiv = ui.div([dfInput, errorDiv]);
  if (Object.keys(dataSourceFunctionsMap).length === 0) {
  // functions that have special tag and are applicable for data source. they should return a dataframe with molecules
    const dataSourceFunctions =
      Array.from(new Set(DG.Func.find({meta: {role: dataSourceTag}}).concat(DG.Func.find({tags: [dataSourceTag]}))));
    // for display purposes we use friendly name of the function
    dataSourceFunctions.forEach((func) => {
      dataSourceFunctionsMap[func.friendlyName ?? func.name] = func;
    });
  }
  let funcCall: DG.FuncCall | null = null;
  let funcCallScope: Scope | undefined;
  // each data source function can have some parameters like for example number of rows to return
  // whenever user selects a function we create a FuncCall object and get an editor for it
  const onDataFunctionChange = async () => {
    const func = dataSourceFunctionsMap[dataSourceFunctionInput.value!];
    funcCall = func.prepare();
    funcCallScope?.dispose();
    funcCallScope = options.owner.child('data-source-call');
    for (const param of funcCall.inputParams.values())
      funcCallScope.subscribe(param.onChanged, edited);
    const editor = await funcCall.getEditor();
    editor.classList.remove('ui-form');
    funcEditorDiv.innerHTML = '';
    funcEditorDiv.appendChild(editor);
    edited();
  };
  const funcEditorDiv = ui.div([]);
  const dataSourceFunctionInput = ui.input.choice(
    C.i18n.dataSourceFunction, {value: template.queryFunctionName ?? Object.keys(dataSourceFunctionsMap)[0],
      items: Object.keys(dataSourceFunctionsMap), onValueChanged: onDataFunctionChange});
  // call the onchange function to create an editor for the first function
  await onDataFunctionChange();
  if (template.queryFunctionName)
    dataSourceFunctionInput.root.getElementsByTagName('select').item(0)?.setAttribute('disabled', 'true');
  const functionInputDiv = ui.div([dataSourceFunctionInput, funcEditorDiv]);
  // if the file source is selected as 'File', no other inputs are needed so we hide the function editor
  functionInputDiv.style.display = 'none';
  const dataInputsDiv = ui.div([fileInputDiv, functionInputDiv]);

  const campaignNameInput = ui.input.string('Campaign Name', {value: '', nullable: true});
  campaignNameInput.setTooltip('Optional friendly name for the campaign');

  // campaign properties. each template might have number of additional fields that should
  // be filled by user for the campaign. they are cast into DG.Property objects and displayed as a form
  const campaignProps = template.campaignFields
    .map((field) => field.type === DG.SEMTYPE.MOLECULE ?
      ({...field, type: 'String', semtype: DG.SEMTYPE.MOLECULE}) : field)
    .map((field) =>
      DG.Property.fromOptions(
        {name: field.name, type: CampaignFieldTypes[field.type as keyof typeof CampaignFieldTypes],
          nullable: !field.required, ...(field.semtype ? {semType: field.semtype} : {})}));
  const campaignPropsObject: {[key: string]: any} = {};
  const campaignPropsForm = campaignProps.length ?
    ui.input.form(campaignPropsObject, campaignProps, {onValueChanged: edited}) : ui.div();
  campaignPropsForm.classList.remove('ui-form');
  // displaying function editor or file input depending on the data source type
  if (template.dataSourceType === 'File') {
    fileInputDiv.style.display = 'inherit';
    functionInputDiv.style.display = 'none';
  } else {
    fileInputDiv.style.display = 'none';
    functionInputDiv.style.display = 'inherit';
  }

  const isEmpty = (value: unknown) => value === null || value === undefined || value === '';
  const problem = computed(() => {
    edits.value;
    if (template.dataSourceType === 'File' && !fileDf)
      return 'Pick a file';
    if (template.dataSourceType !== 'File') {
      const param = [...funcCall?.inputParams.values() ?? []].find((p) => !p.property.nullable && isEmpty(p.value));
      if (param)
        return `${param.property.caption || param.name}: ${EMPTY}`;
    }
    const field = template.campaignFields.find((f) => f.required && isEmpty(campaignPropsObject[f.name]));
    return field ? `${field.name}: ${EMPTY}` : null;
  });

  const dataOf = async (): Promise<INewCampaignResult> => {
    const friendlyName = campaignNameInput.value?.trim() || undefined;
    const extraData = options.getExtraData?.();
    if (template.dataSourceType === 'File')
      return {df: fileDf!, type: 'File', campaignProps: campaignPropsObject, friendlyName, extraData};
    const funcCallInputs: {[key: string]: any} = {};
    for (const [key, value] of Object.entries(funcCall!.inputs))
      funcCallInputs[key] = value;
    const df: DG.DataFrame = await dataSourceFunctionsMap[dataSourceFunctionInput.value!].apply(funcCallInputs);
    return {df, type: 'Query', campaignProps: campaignPropsObject, friendlyName, extraData};
  };
  const busy = signal(false);
  const startCampaignButton = options.owner.build(() => button(C.i18n.StartCampaign, async () => {
    busy.value = true;
    try {
      const data = await dataOf();
      if (data.df.rowCount === 0) {
        grok.shell.warning(`${template.dataSourceType === 'File' ? 'The file has' : `${dataSourceFunctionInput.value} returned`} ` +
          'no rows, so no campaign was created');
        return;
      }
      await options.start(data);
    } finally {
      busy.value = false;
    }
  }, {primary: true, disabled: computed(() => busy.value ? 'Starting…' : problem.value)}));

  const form = ui.div([
    campaignNameInput.root,
    ...(options.extraFormElements ?? []),
    dataInputsDiv,
    ...(campaignProps.length ? [campaignPropsForm] : [])]);
  const buttons = ui.buttonsInput([]);
  buttons.appendChild(startCampaignButton);
  form.appendChild(buttons);
  // the platform form around the rows, as the landing had it: one label column, START at the form's end
  const root = ui.div([form], 'ui-form');
  const column = options.column;
  if (column && typeof ResizeObserver !== 'undefined') {
    const observer = new ResizeObserver(() => column.report(root, Math.max(0, ...Array.from(root
      .querySelectorAll<HTMLElement>('label.ui-input-label')).filter((label) => label.textContent).map((label) => label.offsetWidth))));
    observer.observe(root);
    options.owner.own(() => {
      observer.disconnect();
      column.leave(root);
    });
  }
  return root;
}
