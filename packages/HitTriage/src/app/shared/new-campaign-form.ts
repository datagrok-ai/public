/* The New Campaign form of Hit Design and PeptiHit: Campaign Name and an input per field of the template, on the
   picker's label column, then START — disabled while a required field is empty, with the reason as its tooltip and
   in the line beside it, which otherwise says what START will make. */
import {AsyncValue, BoolInput, Control, DateInput, Form, Input, NumberInput, TextInput, button, component, computed,
  div, signal, statusLine} from '@datagrok-libraries/u2';
import type {LabelColumn} from '@datagrok-libraries/u2';
import {moleculeInput} from '@datagrok-libraries/u2/dg';
import {HitDesignTemplate, HitTriageCampaignField} from '../types';
import {i18n} from '../consts';
import {PICKER_PX} from './template-picker';

export interface NewCampaignRequest {
  /** The field values as campaign.json keeps them: an empty field is left out, a date as {_value, isUtc, date}. */
  fields: {[name: string]: unknown};
  name?: string;
}

export interface NewCampaignFormOptions {
  template: HitDesignTemplate;
  column: LabelColumn;
  /** The code the campaign is likely to get ("KHL-13"); the save takes the one that is free then. */
  nextCode: () => Promise<string>;
  start: (request: NewCampaignRequest) => Promise<unknown>;
}

type FieldInput = readonly [HitTriageCampaignField, Input<any>];

export interface NewCampaignForm {
  readonly form: Form;
}

export function newCampaignForm(options: NewCampaignFormOptions): Control & NewCampaignForm {
  return component(() => {
    const form = new Form({labelColumn: options.column, maxWidth: PICKER_PX});
    const name = new TextInput({label: 'Campaign Name', name: 'campaign-name', placeholder: 'Optional',
      tooltipText: 'New campaign name. If empty, campaign code will be used.'});
    form.add(name);
    const fields: FieldInput[] = options.template.campaignFields.map((f) => [f, fieldInput(f)] as const);
    for (const [, input] of fields)
      form.add(input);
    const code = new AsyncValue(() => options.nextCode(), {debounceMs: 0});
    const busy = signal(false);
    const missing = computed(() => fields.filter(([f, input]) => f.required && input.validity.value !== null)
      .map(([f]) => f.name));
    const start = button(i18n.StartCampaign, async () => {
      busy.value = true;
      try {
        await options.start({fields: fieldValues(fields), name: name.peek()?.trim() || undefined});
      } finally {
        busy.value = false;
      }
    }, {primary: true, disabled: computed(() => busy.value ? 'Starting…' : form.problem.value)});
    const quiet = computed(() => missing.value.length > 0 ? `Fill in ${missing.value.join(', ')} to start` :
      code.value === undefined ? '' : `Creates ${code.value} from ${options.template.name}`);
    // in red once a field shows its error, named by its caption ("Owner: Value can't be empty")
    form.addButtons((row) => row.append(start, statusLine(form.shownProblem, quiet)), {align: 'editors'});
    return {root: div([form.root]), form};
  }, {name: 'new-campaign'});
}

/** The editor of a template field: text, number, date, checkbox, or the sketcher for a molecule. */
function fieldInput(field: HitTriageCampaignField): Input<any> {
  const options = {label: field.name, name: field.name, required: field.required === true};
  if (field.type === 'Molecule' || field.semtype === 'Molecule')
    return moleculeInput(field.name, {name: field.name, required: field.required === true});
  switch (field.type) {
  case 'Number':
    return new NumberInput(options);
  case 'Date':
    return new DateInput(options);
  case 'Boolean':
    return new BoolInput({label: field.name, name: field.name});
  default:
    return new TextInput(options);
  }
}

/** What the inputs hold, as campaign.json keeps it: a Date as {_value, isUtc, date}, the fields left empty (a
 * Boolean left unchecked) out, as the platform form left them. */
function fieldValues(fields: readonly FieldInput[]): {[name: string]: unknown} {
  const values: {[name: string]: unknown} = {};
  for (const [field, input] of fields) {
    const value = input.peek();
    if (Input.isEmpty(value) || value === false)
      continue;
    values[field.name] = value instanceof Date ?
      {_value: value.getTime(), isUtc: false, date: value.toISOString()} : value;
  }
  return values;
}
