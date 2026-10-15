/* The kind catalog: one entry per `data-u2` value u2 stamps — what a test calls it, where its
   name and parts are, which state keys it fills and what a user does to it. Data only: the manifest
   generator merges it with the registry into the manifest's `kinds`, and the bdd library generates its
   vocabulary from that. Contract: core/docs/features/ui2/overhaul/CONTRACT.md §4–§5. A new `data-u2`
   value gets its entry here (tests/kinds.test.js fails until it has one). */

export type U2Family = 'input' | 'action' | 'container' | 'collection' | 'navigation' | 'display' |
  'overlay' | 'composite';

/** What the user does to a control or an item — the bdd library's gesture verbs. */
export type U2Gesture =
  | 'type' | 'clear' | 'pick' | 'toggle' | 'press' | 'expand' | 'collapse' | 'select'
  | 'open' | 'close' | 'drag' | 'resize' | 'hover' | 'scroll' | 'add' | 'remove' | 'rename';

/** The keys of `IControlState` (js-api `u2core/protocol.ts`). */
export type U2StateKey = 'kind' | 'tag' | 'name' | 'label' | 'enabled' | 'visible' | 'focused' |
  'readOnly' | 'required' | 'valid' | 'error' | 'busy' | 'value' | 'text' | 'placeholder' | 'open' |
  'expanded' | 'checked' | 'pressed' | 'selected' | 'items' | 'selection' | 'current' | 'extra';

export type U2NameSource = 'label' | 'aria' | 'text' | 'title' | 'placeholder';

export interface U2KindEntry {
  /** The `data-u2` value. Unique. */
  kind: string;
  /** The bdd noun it answers to ('text input'); several kinds may share a noun (number/bigint/qnum). */
  noun: string;
  aliases?: string[];
  family: U2Family;
  /** Registry tag building it, when one does. */
  tag?: string;
  /** ARIA role of the root. */
  role?: string;
  /** ARIA role of the focus target (inputs: 'textbox', 'combobox', 'checkbox', 'slider' …). */
  focusRole?: string;
  /** Where the user-visible name comes from, in priority order. */
  nameFrom: U2NameSource[];
  /** Part name → selector relative to the root. Always includes the shared input parts for inputs. */
  parts: Record<string, string>;
  /** State keys this kind fills (beyond the always-present kind/enabled/visible/focused). */
  state: U2StateKey[];
  /** Settable props: the tag's registry prop names minus appearance and `structural` ones. Empty for
   * tagless kinds. */
  props: string[];
  /** Present on collection-like kinds. */
  items?: {role: string; noun: string; gestures: U2Gesture[]; multi?: boolean};
  /** Gestures on the control itself. */
  gestures: U2Gesture[];
  /** Legacy sub-element ids stamped as `data-u2` (fb-search, ff-run …): the kind they are a part of.
   * Such an entry's `noun` is its part name, not a vocabulary noun. */
  partOf?: string;
}

/** An entry as written here; the generator adds `tag` and `props` from the registry. */
export type U2KindDef = Omit<U2KindEntry, 'tag' | 'props'>;

/** The manifest's `kinds`: every definition with the registry tag building it and that tag's
 * settable props, keys in one order; a `partOf` entry also becomes a part of its owner under its
 * noun. Shared by the u2 and ems generators and their drift tests. */
export function kindCatalog(defs: readonly U2KindDef[], tagOf: (kind: string) => string | undefined,
  propsOf: (tag: string) => string[]): U2KindEntry[] {
  const children = new Map<string, Record<string, string>>();
  for (const def of defs) {
    if (def.partOf !== undefined)
      children.set(def.partOf, {...children.get(def.partOf), [def.noun]: `[data-u2="${def.kind}"]`});
  }
  return defs.map((def) => {
    const tag = tagOf(def.kind);
    const entry: U2KindEntry = {kind: def.kind, noun: def.noun, aliases: def.aliases, family: def.family, tag,
      role: def.role, focusRole: def.focusRole, nameFrom: def.nameFrom,
      parts: {...def.parts, ...children.get(def.kind)}, state: def.state,
      props: tag === undefined ? [] : propsOf(tag), items: def.items, gestures: def.gestures, partOf: def.partOf};
    for (const key of Object.keys(entry) as (keyof U2KindEntry)[]) {
      if (entry[key] === undefined)
        delete entry[key];
    }
    return entry;
  });
}

const INPUT_PARTS: Record<string, string> = {
  label: '[data-u2-part="label"]',
  editor: '[data-u2-part="editor"]',
  options: '[data-u2-part="options"]',
  error: '[data-u2-part="error"]',
};

const NAME_FROM: Record<U2Family, U2NameSource[]> = {
  input: ['label', 'aria', 'placeholder'],
  action: ['text', 'aria', 'title'],
  container: ['title', 'aria'],
  collection: ['aria'],
  navigation: ['aria'],
  display: ['text', 'aria'],
  overlay: ['title', 'aria'],
  composite: ['aria'],
};

type Extra = Partial<Omit<U2KindDef, 'kind' | 'noun' | 'family'>>;

function def(family: U2Family, kind: string, noun: string, extra: Extra = {}): U2KindDef {
  return {kind, noun, family, nameFrom: NAME_FROM[family], parts: {}, state: [], gestures: [], ...extra};
}

/** An input: the shared parts, and its label is always read. */
function input(kind: string, noun: string, focusRole: string | undefined, state: U2StateKey[],
  gestures: U2Gesture[], extra: Extra = {}): U2KindDef {
  return def('input', kind, noun, {focusRole, parts: INPUT_PARTS, state: ['label', ...state], gestures, ...extra});
}

function part(kind: string, name: string, owner: string, family: U2Family): U2KindDef {
  return {kind, noun: name, family, nameFrom: [], parts: {}, state: [], gestures: [], partOf: owner};
}

const ITEM_PICK: NonNullable<U2KindDef['items']> = {role: 'option', noun: 'item', gestures: ['pick']};

const TEXT_STATE: U2StateKey[] = ['value', 'text', 'placeholder', 'valid', 'error', 'readOnly', 'required'];
const NUMBER_STATE: U2StateKey[] = ['value', 'text', 'valid', 'error', 'readOnly', 'required'];
/** What every `Input` reports (input-base.ts getState). */
const INPUT_STATE: U2StateKey[] = ['value', 'text', 'valid', 'error', 'readOnly', 'required'];
/** An `Input` with a popup of its own. */
const PICKER_STATE: U2StateKey[] = [...INPUT_STATE, 'open'];
const COMBO_STATE: U2StateKey[] = ['value', 'text', 'open', 'placeholder'];

export const U2_KINDS: readonly U2KindDef[] = [
  // --- inputs ---------------------------------------------------------------------------------------
  input('text-input', 'text input', 'textbox', TEXT_STATE, ['type', 'clear'], {aliases: ['text field', 'textbox']}),
  input('text-area', 'text area', 'textbox', TEXT_STATE, ['type', 'clear'], {aliases: ['textarea', 'multiline input']}),
  input('number-input', 'number input', 'spinbutton', NUMBER_STATE, ['type', 'clear'],
    {aliases: ['numeric input', 'number field']}),
  input('bigint-input', 'number input', 'spinbutton', NUMBER_STATE, ['type', 'clear']),
  input('qnum-input', 'number input', 'textbox', NUMBER_STATE, ['type', 'clear']),
  input('slider-input', 'slider', 'slider', INPUT_STATE, ['drag'], {aliases: ['slider input']}),
  input('range-slider', 'range slider', 'slider', ['value', 'text'], ['drag'], {aliases: ['range input']}),
  input('bool-input', 'checkbox', 'checkbox', [...INPUT_STATE, 'checked'], ['toggle'],
    {aliases: ['bool input', 'switch', 'toggle']}),
  input('choice-input', 'choice input', 'combobox', [...INPUT_STATE, 'items', 'selection'],
    ['open', 'pick'], {aliases: ['dropdown', 'choice', 'select'], items: ITEM_PICK}),
  input('radio-input', 'radio input', 'radiogroup', [...INPUT_STATE, 'items'], ['pick'],
    {aliases: ['radio group', 'radio'], items: {role: 'radio', noun: 'radio', gestures: ['pick']}}),
  input('multi-choice-input', 'multi choice input', 'group', [...INPUT_STATE, 'items', 'selection'], ['toggle'],
    {aliases: ['multi choice'], items: {role: 'checkbox', noun: 'checkbox', gestures: ['toggle'], multi: true}}),
  input('chips-input', 'chips input', 'group', [...INPUT_STATE, 'items', 'selection'], ['toggle'],
    {items: {role: 'button', noun: 'chip', gestures: ['toggle'], multi: true}}),
  input('tags-input', 'tags input', 'textbox', [...PICKER_STATE, 'items'], ['type', 'add', 'remove'],
    {aliases: ['tags field'], items: {role: 'listitem', noun: 'tag', gestures: ['remove']}}),
  input('list-input', 'list input', 'textbox', [...INPUT_STATE, 'items'], ['type', 'add', 'remove'],
    {items: {role: 'listitem', noun: 'tag', gestures: ['remove']}}),
  input('map-input', 'map input', 'grid', [...INPUT_STATE, 'items'], ['add', 'remove', 'type'],
    {aliases: ['key value input'], items: {role: 'row', noun: 'item', gestures: ['remove']}}),
  input('color-input', 'color input', 'textbox', INPUT_STATE, ['type', 'pick'], {aliases: ['color picker']}),
  input('icon-input', 'icon input', 'button', PICKER_STATE, ['open', 'pick'], {aliases: ['icon picker']}),
  input('font-input', 'font input', 'textbox', [...PICKER_STATE, 'items', 'selection'], ['open', 'pick'],
    {aliases: ['font picker'], items: {role: 'option', noun: 'size', gestures: ['pick']}}),
  input('image-input', 'image input', undefined, INPUT_STATE, ['clear'], {aliases: ['image picker']}),
  input('date-input', 'date input', 'textbox', PICKER_STATE, ['type', 'open', 'pick'],
    {aliases: ['date field', 'datetime input', 'date picker']}),
  input('datetime-input', 'date input', 'textbox', PICKER_STATE, ['type', 'open', 'pick']),
  input('combobox', 'combobox', 'combobox', [...COMBO_STATE, 'busy', 'items'], ['type', 'open', 'pick'],
    {aliases: ['combo box', 'editable choice'], items: ITEM_PICK}),
  input('typeahead', 'suggest input', 'combobox', [...COMBO_STATE, 'busy', 'items'], ['type', 'open', 'pick'],
    {aliases: ['autocomplete', 'typeahead'], items: ITEM_PICK}),
  input('suggest-input', 'suggest input', 'combobox', [...PICKER_STATE, 'placeholder', 'busy', 'items'],
    ['type', 'open', 'pick'], {items: ITEM_PICK}),
  input('multi-select', 'multi select', 'combobox', [...PICKER_STATE, 'items', 'selection'], ['open', 'toggle'],
    {aliases: ['multi select input'], items: {role: 'option', noun: 'item', gestures: ['toggle'], multi: true}}),
  input('access-grid', 'access grid', 'grid', [...INPUT_STATE, 'items'], ['add', 'toggle', 'remove'],
    {items: {role: 'row', noun: 'table row', gestures: ['toggle', 'remove']}}),
  input('message-input', 'message input', 'textbox', ['value', 'text', 'open', 'placeholder'], ['type'],
    {aliases: ['prompt input', 'chat input']}),
  input('filter-builder', 'filter builder', 'group', ['value', 'valid', 'error'], ['add', 'remove', 'type'], {
    parts: {...INPUT_PARTS, 'mode': '[data-u2-part="mode"]', 'none': '[data-u2-part="none"]',
      'flatten': '[data-u2-part="flatten"]', 'hint mode': '[data-u2-part="hint-mode"]',
      'hint': '[data-u2-part="hint"]', 'query': '[data-u2-part="query"]'},
    items: {role: 'row', noun: 'filter row', gestures: ['remove']},
  }),
  input('filter-query-input', 'filter query', 'textbox', ['value', 'valid', 'error'], ['type'],
    {aliases: ['filter query input']}),
  input('ref-input', 'ref input', 'combobox', COMBO_STATE, ['type', 'open', 'pick'], {items: ITEM_PICK}),
  input('column-combo', 'column input', 'combobox', [...PICKER_STATE, 'items', 'selection'], ['open', 'pick'],
    {aliases: ['column picker', 'column selector', 'column combobox'], items: ITEM_PICK}),
  input('column-picker', 'column input', 'combobox', [...PICKER_STATE, 'items'], ['open', 'pick'], {items: ITEM_PICK}),
  input('columns-input', 'columns input', 'button', [...PICKER_STATE, 'items', 'selection'], ['open', 'toggle'],
    {aliases: ['columns picker'], items: {role: 'option', noun: 'item', gestures: ['toggle'], multi: true}}),
  input('columns-map-input', 'columns input', 'group', INPUT_STATE, ['open', 'toggle']),
  input('aggregated-columns-input', 'columns input', 'group', INPUT_STATE, ['open', 'toggle']),
  input('file-input', 'file input', 'group', [...INPUT_STATE, 'busy'], ['pick'],
    {aliases: ['file field', 'files input']}),
  input('files-input', 'file input', 'group', INPUT_STATE, ['pick']),
  input('function-input', 'function input', 'button', [...PICKER_STATE, 'placeholder'], ['open', 'pick'],
    {aliases: ['func input']}),
  input('func-call-input', 'function input', 'button', PICKER_STATE, ['open', 'pick']),
  input('rsa-input', 'rsa input', 'textbox', INPUT_STATE, ['type'], {aliases: ['key input']}),
  input('dynamic-input', 'dynamic input', undefined, INPUT_STATE, ['type', 'pick']),
  input('dart-input', 'dart input', undefined, INPUT_STATE, ['type'], {aliases: ['bridged input']}),
  def('input', 'readonly-field', 'read-only field', {nameFrom: ['label'],
    parts: {label: '.u2-input-label', value: '[data-u2-part="readonly-value"]'},
    state: ['label', 'value', 'text', 'readOnly']}),

  // --- actions --------------------------------------------------------------------------------------
  def('action', 'button', 'button', {role: 'button', state: ['label', 'busy', 'pressed'], gestures: ['press']}),
  def('action', 'icon-button', 'icon button', {role: 'button', nameFrom: ['aria', 'title'],
    state: ['label', 'pressed'], gestures: ['press']}),
  def('action', 'dropdown-button', 'dropdown button', {aliases: ['menu button', 'split button'], role: 'button',
    state: ['label', 'open'], gestures: ['press', 'open', 'pick']}),
  def('action', 'button-group', 'button group', {aliases: ['segmented control'], role: 'group', nameFrom: ['aria'],
    state: ['value', 'current', 'items'], items: {role: 'button', noun: 'button', gestures: ['press']},
    gestures: ['pick']}),
  def('action', 'row-actions', 'row actions', {aliases: ['actions'], role: 'toolbar', nameFrom: ['aria'],
    state: ['items'], items: {role: 'button', noun: 'button', gestures: ['press']}, gestures: ['press']}),

  // --- containers -----------------------------------------------------------------------------------
  // a container's parts are anchored at its root (`:scope > …`), so a section or a card it holds
  // does not answer for its title, body or footer
  def('container', 'section', 'section', {aliases: ['pane', 'accordion pane'], role: 'group',
    parts: {title: ':scope > [data-u2-part="header"] > [data-u2-part="title"]',
      header: ':scope > [data-u2-part="header"]', body: ':scope > [data-u2-part="body"]'},
    state: ['label', 'expanded'], gestures: ['expand', 'collapse']}),
  def('container', 'accordion', 'accordion', {nameFrom: ['aria'], state: ['items'],
    items: {role: 'button', noun: 'accordion header', gestures: ['expand', 'collapse']},
    gestures: ['expand', 'collapse']}),
  def('container', 'tabs', 'tabs', {aliases: ['tab strip', 'tab control'], nameFrom: ['aria'],
    parts: {'close button': '[data-u2-part="close"]'}, state: ['current', 'items', 'selection'],
    items: {role: 'tab', noun: 'tab', gestures: ['select']}, gestures: ['select']}),
  def('container', 'wizard', 'wizard', {
    nameFrom: ['aria'],
    parts: {'steps': ':scope > [data-u2-part="steps"]', 'body': ':scope > [data-u2-part="body"]',
      'footer': ':scope > [data-u2-part="footer"]',
      'status': ':scope > [data-u2-part="footer"] > * > [data-u2-part="status"]',
      'reason': ':scope > [data-u2-part="footer"] > * > [data-u2-part="reason"]',
      'error row': ':scope > [data-u2-part="footer"] > * > [data-u2-part="error-row"]'},
    state: ['current', 'busy', 'valid', 'error', 'items'],
    items: {role: 'listitem', noun: 'wizard step', gestures: ['select']}, gestures: ['select'],
  }),
  def('container', 'splitter', 'splitter', {aliases: ['split panel'], nameFrom: ['aria'],
    parts: {sash: '[data-u2-part="sash"]'}, state: ['extra'], gestures: ['drag', 'resize']}),
  def('container', 'card', 'card', {
    nameFrom: ['title', 'text'],
    parts: {title: ':scope > [data-u2-part="header"] > * > [data-u2-part="title"]',
      subtitle: ':scope > [data-u2-part="header"] > * > [data-u2-part="subtitle"]',
      header: ':scope > [data-u2-part="header"]', media: ':scope > [data-u2-part="media"]',
      body: ':scope > [data-u2-part="body"]', footer: ':scope > [data-u2-part="footer"]',
      actions: ':scope > [data-u2-part="header"] > [data-u2-part="actions"]'},
    state: ['label', 'selected'], gestures: ['press'],
  }),
  def('container', 'form', 'form', {role: 'form', nameFrom: ['aria'], state: ['valid', 'error', 'items']}),
  def('container', 'object-form', 'object form', {role: 'form', nameFrom: ['aria'],
    state: ['valid', 'error', 'items']}),
  def('container', 'func-form', 'function form', {aliases: ['func form'], role: 'form', nameFrom: ['aria'],
    state: ['valid', 'error', 'busy', 'items']}),
  // layout only: never addressed by what they show, only named (CONTRACT §5.3)
  def('container', 'panel', 'layout panel', {nameFrom: ['aria']}),
  def('container', 'div-v', 'vertical stack', {nameFrom: ['aria']}),
  def('container', 'div-h', 'horizontal stack', {nameFrom: ['aria']}),
  def('container', 'async-view', 'async view', {aliases: ['loading view'], nameFrom: ['aria'],
    state: ['busy', 'error']}),

  // --- collections ----------------------------------------------------------------------------------
  def('collection', 'list', 'list', {role: 'listbox', state: ['label', 'items', 'selection'],
    items: {role: 'option', noun: 'item', gestures: ['select']}, gestures: ['select', 'scroll']}),
  def('collection', 'tree', 'tree', {role: 'tree', state: ['items', 'selection'],
    parts: {'item-label': '[data-u2-part="item-label"]'},
    items: {role: 'treeitem', noun: 'tree node', gestures: ['select', 'expand', 'collapse', 'toggle']},
    gestures: ['select', 'expand', 'collapse']}),
  def('collection', 'table', 'table', {role: 'table', state: ['items', 'selection'],
    parts: {'item-label': '[data-u2-part="item-label"]', 'section': '[data-u2-part="section"]'},
    items: {role: 'row', noun: 'table row', gestures: ['select']}, gestures: ['select', 'scroll']}),
  def('collection', 'data-table', 'table', {role: 'grid', state: ['items', 'selection'],
    parts: {'item-label': '[data-u2-part="item-label"]'},
    items: {role: 'row', noun: 'table row', gestures: ['select']}, gestures: ['select', 'scroll']}),
  def('collection', 'grid', 'grid', {aliases: ['virtual grid', 'cell grid'], role: 'listbox',
    state: ['items', 'selection'], items: {role: 'option', noun: 'cell', gestures: ['select', 'pick']},
    gestures: ['select', 'scroll']}),
  def('collection', 'property-grid', 'property grid', {aliases: ['property editor'], role: 'grid', state: ['items'],
    parts: {'item-label': '[data-u2-part="item-label"]', 'category': '[data-u2-part="category"]'},
    items: {role: 'row', noun: 'property', gestures: ['type']}, gestures: ['type']}),
  // two forms (identity, typed fields): its items are the fields, found by their labels like a form's
  def('collection', 'property-editor', 'property grid', {state: ['items'], gestures: ['type']}),
  def('collection', 'functions-browser', 'functions browser', {state: ['items', 'selection'],
    items: {role: 'option', noun: 'item', gestures: ['select']}, gestures: ['type', 'select']}),
  def('collection', 'func-call-history-browser', 'history browser',
    {aliases: ['func call history browser', 'run history'], state: ['items', 'busy', 'error'],
      items: {role: 'option', noun: 'item', gestures: ['select']}, gestures: ['select']}),

  // --- navigation -----------------------------------------------------------------------------------
  def('navigation', 'breadcrumbs', 'breadcrumbs', {aliases: ['breadcrumb bar'], role: 'navigation',
    state: ['current', 'items'], items: {role: 'listitem', noun: 'breadcrumb', gestures: ['press']},
    gestures: ['press']}),
  def('navigation', 'menu', 'menu', {role: 'menu', state: ['open', 'items'],
    items: {role: 'menuitem', noun: 'menu item', gestures: ['pick', 'hover']}, gestures: ['pick', 'hover']}),
  def('navigation', 'menu-bar', 'menu bar', {role: 'menubar', state: ['items'],
    items: {role: 'menuitem', noun: 'menu item', gestures: ['open', 'pick']}, gestures: ['open', 'pick']}),
  def('navigation', 'toolbar', 'toolbar', {role: 'toolbar', state: ['items'],
    items: {role: 'button', noun: 'button', gestures: ['press']}, gestures: ['press']}),

  // --- overlays -------------------------------------------------------------------------------------
  def('overlay', 'dialog', 'dialog', {role: 'dialog',
    parts: {'title': ':scope > .u2-dialog-title > [data-u2-part="title"]',
      'close button': ':scope > .u2-dialog-title > [data-u2-part="close"]',
      'body': ':scope > [data-u2-part="body"]', 'footer': ':scope > [data-u2-part="footer"]',
      'error row': ':scope > [data-u2-part="error-row"]',
      'spinner': ':scope > [data-u2-part="footer"] > * > [data-u2-part="busy"]'},
    state: ['label', 'open', 'busy', 'error', 'valid'], gestures: ['close']}),
  def('overlay', 'tooltip', 'tooltip', {role: 'tooltip', nameFrom: ['text'], state: ['text'], gestures: ['hover']}),
  def('overlay', 'notify', 'notification', {aliases: ['balloon', 'toast', 'message'], role: 'status',
    nameFrom: ['text'], state: ['text', 'extra'], gestures: ['close']}),
  def('overlay', 'tour', 'tour', {role: 'dialog', nameFrom: ['title'], state: ['current', 'items'],
    gestures: ['press']}),

  // --- display --------------------------------------------------------------------------------------
  def('display', 'icon', 'icon', {nameFrom: ['aria'], state: ['label', 'extra']}),
  def('display', 'badge', 'badge', {state: ['text']}),
  def('display', 'count-badge', 'badge', {state: ['text']}),
  def('display', 'dot', 'badge', {state: ['text']}),
  def('display', 'tag', 'tag', {nameFrom: ['text'], state: ['text', 'selected', 'pressed']}),
  def('display', 'chip', 'chip', {nameFrom: ['text'], state: ['text', 'selected', 'pressed']}),
  def('display', 'bind-chip', 'chip', {nameFrom: ['text'], state: ['text', 'selected', 'pressed']}),
  def('display', 'msg-chip', 'chip', {nameFrom: ['text'], state: ['text', 'selected', 'pressed']}),
  def('display', 'stat-card', 'stat card', {aliases: ['stat'], nameFrom: ['label', 'text'],
    parts: {label: '[data-u2-part="label"]', value: '[data-u2-part="value"]', delta: '[data-u2-part="delta"]'},
    state: ['label', 'value', 'text', 'extra']}),
  def('display', 'progress-bar', 'progress bar', {aliases: ['progress'], role: 'progressbar',
    nameFrom: ['label', 'aria'], parts: {label: '[data-u2-part="label"]', track: '[data-u2-part="track"]',
      percent: '[data-u2-part="percent"]'}, state: ['value', 'text', 'busy']}),
  def('display', 'entity-card', 'entity card', {aliases: ['entity'], nameFrom: ['title', 'text'],
    parts: {title: '[data-u2-part="title"]'}, state: ['label']}),

  // --- composite (designer, filter builder internals) ------------------------------------------------
  def('composite', 'designer', 'designer', {state: ['busy']}),
  def('composite', 'palette', 'palette'),
  def('composite', 'tray', 'tray'),
  def('composite', 'filter-group', 'filter group', {parts: {'connector': '[data-u2-part="connector"]',
    'not': '[data-u2-part="not"]', 'add': '[data-u2-part="add"]', 'add group': '[data-u2-part="add-group"]',
    'rows': '[data-u2-part="rows"]', 'handle': '[data-u2-part="handle"]'}, gestures: ['add', 'remove', 'drag']}),
  def('composite', 'filter-row', 'filter row', {parts: {property: '[data-u2-part="prop"]',
    operator: '[data-u2-part="op"]', value: '[data-u2-part="value"]', remove: '[data-u2-part="remove"]',
    handle: '[data-u2-part="handle"]'}, gestures: ['type', 'pick', 'remove', 'drag']}),

  // --- legacy sub-element ids: parts of their owner -------------------------------------------------
  part('fb-search', 'search', 'functions-browser', 'collection'),
  part('fb-list', 'list', 'functions-browser', 'collection'),
  part('fb-status', 'status', 'functions-browser', 'collection'),
  part('fb-clear', 'clear button', 'functions-browser', 'collection'),
  part('fb-roles', 'roles', 'functions-browser', 'collection'),
  part('fb-tags', 'tags', 'functions-browser', 'collection'),
  part('fb-empty', 'empty message', 'functions-browser', 'collection'),
  part('fb-panes', 'panes', 'functions-browser', 'collection'),
  part('fb-role', 'role', 'functions-browser', 'collection'),
  part('fb-tag', 'tag', 'functions-browser', 'collection'),
  part('fch-list', 'list', 'func-call-history-browser', 'collection'),
  part('fch-retry', 'retry button', 'func-call-history-browser', 'collection'),
  part('fch-state', 'state', 'func-call-history-browser', 'collection'),
  part('ff-run', 'run button', 'func-form', 'container'),
  part('ff-history-icon', 'history icon', 'func-form', 'container'),
];
