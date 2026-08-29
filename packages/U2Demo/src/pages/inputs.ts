import {
  signal, computed, Scope,
  divV, divH, span, h3,
  TextInput, TextArea, BoolInput, NumberInput, ChoiceInput, MultiChoiceInput, IconInput,
  SuggestInput, DynamicInput, Combobox, AsyncView, RangeSlider, MultiSelect, ButtonGroup,
} from '@datagrok-libraries/u2';
import {functionInput} from '@datagrok-libraries/u2/src/dg/index.js';
import {readout, CITIES} from './common';
import '../../css/u2demo-inputs.css';

export function basicInputsPage(): HTMLElement {
  const query = signal('');
  const notify = signal(true);

  const name = new TextInput({label: 'Name', value: 'Aspirin', tooltipText: 'Compound name'});
  const search = new TextInput({label: 'Search', search: true, bind: query, placeholder: 'Filter compounds…'});
  const notes = new TextArea({label: 'Notes', placeholder: 'Grows as you type…', autoGrow: true});
  const enabled = new BoolInput({label: 'Enabled', value: true});
  const notifications = new BoolInput({label: 'Notifications', switch: true, bind: notify});
  const dose = new NumberInput({label: 'Dose, mg', mode: 'float', min: 0, max: 1000, step: 0.5, value: 250});
  const replicates = new NumberInput({label: 'Replicates', mode: 'int', min: 1, max: 10, value: 3});
  const stage = new ChoiceInput({label: 'Stage', items: ['Discovery', 'Preclinical', 'Phase I', 'Phase II'],
    value: 'Discovery'});
  const targets = new MultiChoiceInput({label: 'Targets', items: ['GPCR', 'Kinase', 'Ion channel', 'Protease'],
    value: ['Kinase']});

  const glyph = new IconInput({label: 'Icon', value: 'flask'});
  const code = new TextInput({label: 'Code', placeholder: 'required, up to 10 characters'});
  code.addValidator((v) => v.trim().length > 0 ? null : 'Value is required');
  code.addValidator((v) => v.length > 10 ? 'At most 10 characters' : null);

  const scorer = functionInput('Scorer', {placeholder: 'Pick a function…'});
  const city = new SuggestInput({label: 'City', placeholder: 'Type for suggestions…',
    source: async (q) => CITIES.filter((c) => c.toLowerCase().includes(q.toLowerCase()))});
  // DynamicInput shows its value instead of editing it: only code ever writes it, here an effect
  const preview = new DynamicInput<string>({label: 'Preview',
    render: (v) => v ? span(v, 'u2demo-code') : null});
  Scope.ambient!.effect(() => preview.value.value = query.value || null);

  return divV([
    span('Every input is a signal owner: bind to it, or pass your own signal through the bind option.',
      'u2demo-hint'),
    name, search, notes, enabled, notifications, dose, replicates, stage, targets, glyph, code,
    scorer, city, preview,
    readout('search', computed(() => query.value || '(empty)')),
    readout('notifications', computed(() => String(notify.value))),
    readout('dose * replicates', computed(() =>
      String((dose.value.value ?? 0) * (replicates.value.value ?? 0)))),
    readout('code validity', computed(() => code.validity.value ?? 'valid')),
    readout('icon', glyph.value),
    readout('scorer', computed(() => scorer.value.value || '(none)')),
    readout('city', computed(() => city.value.value || '(empty)')),
  ], 'u2demo-page');
}

export function rangeSliderPage(): HTMLElement {
  const lo = signal(20);
  const hi = signal(70);
  const range = new RangeSlider({min: 0, max: 100, step: 1, minRange: 5, lo, hi,
    format: (value) => `${value.toFixed(0)}%`});

  const from = signal(0);
  const to = signal(1000);
  const dose = new RangeSlider({min: 0, max: 1000, step: 10, lo: from, hi: to,
    format: (value) => `${value.toFixed(0)} mg`});

  return divV([
    span('Two handles over one numeric range, each on a signal of its own: bind the pair and the ' +
      'slider is a filter. minRange keeps the handles apart, step quantizes the drag, and format ' +
      'writes the tooltip.', 'u2demo-hint'),
    h3('Percent, minRange 5'), range,
    readout('range', computed(() => `${lo.value.toFixed(0)} … ${hi.value.toFixed(0)}`)),
    h3('Dose, step 10'), dose,
    readout('dose', computed(() => `${from.value.toFixed(0)} … ${to.value.toFixed(0)} mg`)),
  ], 'u2demo-page');
}

export function multiSelectPage(): HTMLElement {
  const targets = new MultiSelect({label: 'Targets', selectAll: true, value: ['Kinase'],
    items: ['GPCR', 'Kinase', 'Ion channel', 'Protease', 'Transporter', 'Nuclear receptor']});

  const actions = new ButtonGroup({items: [
    {id: 'copy', label: 'Copy', icon: 'copy'},
    {id: 'paste', label: 'Paste', icon: 'paste'},
    {id: 'delete', label: 'Delete', icon: 'trash'},
  ]});
  const layout = new ButtonGroup({toggle: 'single', density: 'toolbar', items: [
    {id: 'grid', label: 'Grid'}, {id: 'cards', label: 'Cards'}, {id: 'list', label: 'List'},
  ]});
  layout.selected.value = ['grid'];
  const style = new ButtonGroup({toggle: 'multi', iconOnly: true, density: 'ribbon', items: [
    {id: 'bold', label: 'Bold', icon: 'bold'},
    {id: 'italic', label: 'Italic', icon: 'italic'},
    {id: 'underline', label: 'Underline', icon: 'underline'},
  ]});

  return divV([
    span('Picking several values out of a set: a popup list when the set is long, segmented ' +
      'buttons when it is short enough to show whole. A ButtonGroup with no toggle mode is a ' +
      'plain action bar.', 'u2demo-hint'),
    h3('MultiSelect'), targets,
    readout('targets', computed(() => targets.value.value.join(', ') || '(none)')),
    h3('ButtonGroup'),
    divH([actions, layout, style], 'u2demo-row'),
    readout('layout', computed(() => layout.selected.value.join(', ') || '(none)')),
    readout('style', computed(() => style.selected.value.join(', ') || '(none)')),
  ], 'u2demo-page');
}

export function asyncPage(): HTMLElement {
  const local = new Combobox({items: CITIES, placeholder: 'Pick a city…'});

  const fetchCities = async (q: string, abort: AbortSignal): Promise<string[]> => {
    await new Promise((resolve) => setTimeout(resolve, 400));
    if (abort.aborted)
      return [];
    return CITIES.filter((c) => c.toLowerCase().includes(q.toLowerCase()));
  };

  const remote = new Combobox({items: fetchCities, placeholder: 'Type to search (simulated 400 ms latency)…'});

  const view = AsyncView.owned(fetchCities,
    (items) => divV(items.map((c) => span(c)), 'u2demo-card'),
    {empty: 'No cities match', skeleton: true});
  view.refresh();
  const filter = new TextInput({inline: true, search: true, placeholder: 'Filter the async view…',
    onChanged: (v) => view.refresh(v)});

  return divV([
    span('One async contract — debounce, cancellation, loading/empty/error — shared by every async control.',
      'u2demo-hint'),
    h3('Combobox, local items'), local,
    h3('Combobox, async items'), remote,
    readout('local', local.value), readout('remote', remote.value),
    h3('AsyncView'), filter, view,
  ], 'u2demo-page');
}
