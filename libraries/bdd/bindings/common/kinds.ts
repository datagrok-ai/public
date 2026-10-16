/* Generic kinds: "<qualifier> <kind>" resolves without a registry entry. The u2 nouns come from the
   u2 and EMS manifests (u2-kinds.g.ts, `npm run vocab`): every `data-u2` value the libraries stamp has
   a kind there, and `OVERRIDES` adds only what the manifests cannot know — the Dart client's
   selectors and `name=` conventions (`annotate()` in html_utils.dart), legacy aliases, part names and
   gestures, and the strategy order of a noun whose phrases resolve today — so the same phrase works on
   both UI generations. Inputs share the u2 parts (label, editor, options, error) as `<part> of <input>`.
   The kinds below the u2 nouns are not u2's: the Dart shell, item kinds and generic words. */
import {kind, MatchStrategy} from '../../src/registry.js';
import {registerU2Kinds, type U2Override} from '../../src/u2-vocabulary.js';
import {U2_KINDS} from './u2-kinds.g.js';

const INPUT_LABEL = '[data-u2-part="label"], .u2-input-label, .ui-input-label';
const INPUT_EDITOR = '[data-u2-part="editor"] input, [data-u2-part="editor"] select, [data-u2-part="editor"] textarea, ' +
  '[data-u2-part="editor"] [contenteditable], .ui-input-editor';
const INPUT_PARTS = {
  label: '[data-u2-part="label"], .ui-input-label',
  editor: '[data-u2-part="editor"], .ui-input-editor',
  options: '[data-u2-part="options"]',
  error: '[data-u2-part="error"]',
};
const INPUT_MATCH: MatchStrategy[] = ['name', 'label', 'placeholder', 'dart'];
// the primary text of a u2 row — the name a row goes by when it also shows a signature, a login or a shortcut
const PRIMARY_TEXT = '.u2-fb-label, .u2-tree-label, .u2-typeahead-text, .u2-typeahead-user-name, .u2-multi-select-text, ' +
  '.u2-menu-label, .u2-card-title, .d4-tree-view-node-label';
const ITEM_LABEL = '[data-u2-part="item-label"]';

function u2(...ids: string[]): string {
  return ids.map((id) => `[data-u2="${id}"]`).join(', ');
}

/** A u2 input noun as the Dart client spells it too: its `.ui-input-<type>` root, `input-host-<Caption>`. */
function input(dart: string, aliases: string[] = [], extra: U2Override = {}): U2Override {
  return {...dart ? {selector: dart} : {}, aliases, match: INPUT_MATCH, labelSelector: INPUT_LABEL,
    dartNames: ['input-host-{q}'], editorSelector: INPUT_EDITOR, parts: INPUT_PARTS, ...extra};
}

// a noun whose phrases resolve by data-u2-name today: the manifest's aria source would be tried first
const BY_NAME: U2Override = {match: ['name']};

const OVERRIDES: Record<string, U2Override> = {
  // --- inputs ---
  // a bare <input> outside any input host that names itself (aria-label) is a text input too: the
  // Save project dialog's Name
  'text input': input('.ui-input-text, input[aria-label]:not([type="checkbox"]):not(.ui-input-root *)', ['text field', 'textbox'],
    {match: [...INPUT_MATCH, 'aria']}),
  'text area': input('.ui-input-textarea', ['textarea', 'multiline input']),
  'choice input': input('.ui-input-choice', ['dropdown', 'choice', 'select']),
  'multi choice input': input('', ['multi choice']),
  'number input': input('.ui-input-int, .ui-input-float', ['numeric input', 'number field']),
  // a bare box outside any input host that names itself (aria-label): the Binding site popover's "Show side chains"
  'checkbox': input('.ui-input-bool, .ui-input-bool-switch, .u2-multi-choice-item, .u2-columns-option, ' +
    'input[type="checkbox"][aria-label]:not(.ui-input-root *)', ['bool input', 'switch', 'toggle'],
  {match: [...INPUT_MATCH, 'text', 'aria'], items: false}),
  'date input': input('.ui-input-date', ['date field', 'datetime input', 'date picker']),
  'color input': input('.ui-input-color', ['color picker']),
  'font input': input('', ['font picker']),
  'icon input': input('', ['icon picker']),
  'image input': input('', ['image picker']),
  'list input': input(''),
  'map input': input('', ['key value input']),
  'message input': input('', ['prompt input', 'chat input']),
  'radio input': input('.ui-input-radio', ['radio group', 'radio']),
  'slider': input('.ui-input-slider', ['slider input']),
  'range slider': input('svg[type="range-slider"]', ['range input'], {dartNames: ['input-host-{q}', '{q}']}),
  'suggest input': input('', ['autocomplete', 'typeahead']),
  'combobox': input('.ui-input-editable-choice', ['combo box', 'editable choice']),
  'multi select': input('.ui-input-multi-choice', ['multi select input']),
  'tags input': input('.ui-input-tags', ['tags field']),
  'file input': input('.ui-input-file', ['file field', 'files input']),
  'columns input': input('.ui-input-columns', ['columns picker']),
  'column input': input('.d4-column-selector', ['column picker', 'column selector', 'column combobox'],
    {match: ['name', 'label', 'dart'], dartNames: ['div-column-combobox-{q}', 'div-column-combobox-{q}-', 'input-host-{q}'],
      gestures: {open: 'mousedown'}}),
  'function input': input('', ['func input']),
  'dynamic input': input(''),
  'rsa input': input('', ['key input']),
  'dart input': input('', ['bridged input']),

  // --- actions ---
  'button': {selector: 'button, [role="button"], .u2-btn, .ui-btn', match: ['name', 'text', 'aria', 'dart'],
    dartNames: ['button-{q}'], items: false},
  'icon button': {selector: 'button.u2-icon-btn', dartNames: ['icon-{q}']},
  'row actions': BY_NAME,

  // --- forms and collections ---
  'form': {selector: u2('object-form', 'func-form') + ', .ui-form'},
  'list': {u2Selector: '[data-u2="list"]:not([role="tree"])', selector: '[role="listbox"], [role="list"], .u2-list:not([role="tree"]), .d4-list'},
  'tree': {selector: '[role="tree"], .d4-tree-view'},
  'table': {selector: 'table'},
  'grid': {selector: '.d4-grid'},
  'functions browser': BY_NAME,
  'history browser': BY_NAME,

  // --- display ---
  'icon': {selector: '.u2-icon, .grok-icon, [name^="icon-"], i[class*="fa-"]', match: ['aria', 'name', 'dart'], dartNames: ['icon-{q}', '{q}']},
  'tag': {selector: '.d4-tag', items: false},
  'notification': {u2Selector: '[data-u2="notify"] .u2-notify', selector: '.d4-balloon'},
  'progress bar': {selector: '[role="progressbar"]', labelSelector: '.u2-progress-description'},
  'stat card': {labelSelector: '.u2-stat-label'},
  'tooltip': {selector: '[role="tooltip"], .d4-tooltip'},
  'async view': BY_NAME,
  'tray': BY_NAME,

  // --- navigation ---
  'toolbar': {selector: '[role="toolbar"], .d4-ribbon'},
  // the Dart main menu keeps a permanent [role="menu"] container per group — not a popup
  'menu': {selector: '.d4-menu-popup, [role="menu"]:not(.d4-menu-item-container)'},
  'menu bar': {selector: '[role="menubar"]'},

  // --- containers ---
  // a container's label and parts are its own, never the title, close or footer part of a tab, a
  // section or a card inside it — the manifest's [data-u2-part] selectors reach those too
  'dialog': {selector: '.d4-dialog, [role="dialog"]', match: ['name', 'title', 'dart'],
    labelSelector: '.u2-dialog-title-text, .d4-dialog-title', dartNames: ['dialog-{q}'],
    parts: {title: '.u2-dialog-title, .d4-dialog-title', 'close button': '.u2-dialog-close, .d4-dialog-close',
      footer: '.u2-dialog-buttons, .u2-dialog-footer, .d4-dialog-footer'}},
  'tabs': {selector: '.d4-tab-control'},
  'section': {selector: '.u2-section, .u2-accordion-pane, [name^="div-section--"], .d4-accordion-pane', match: ['name', 'title', 'dart'],
    labelSelector: '.u2-section-title, .u2-section-header, .u2-accordion-title, .d4-accordion-pane-header',
    dartNames: ['div-section--{q}', 'pane-{q}']},
  'accordion': {selector: '.d4-accordion'},
  'card': {selector: u2('stat-card', 'entity-card') + ', .d4-item-card', labelSelector: '.u2-card-title, .card-label'},
  'splitter': BY_NAME,

  // --- entities and the designer ---
  // the chips of a chips input are its items, addressed by key (data-u2-item)
  'chip': {selector: '[data-u2="chips-input"] [data-u2-item]', match: ['text', 'item', 'name']},
  'palette': BY_NAME,
  'designer': BY_NAME,
};

registerU2Kinds(U2_KINDS, OVERRIDES);

// --- inputs ---------------------------------------------------------------------------------------
// a u2 input bridging a Dart one (dart-input) holds the Dart `.ui-input-root`: that inner root is not a
// second input
kind('input', {
  aliases: ['field'],
  selector: '[data-u2$="-input"], ' + u2('combobox', 'multi-select', 'column-picker', 'column-combo', 'range-slider', 'typeahead') +
    ', .ui-input-root:not([data-u2="dart-input"] .ui-input-root)',
  match: INPUT_MATCH,
  labelSelector: INPUT_LABEL,
  dartNames: ['input-host-{q}'],
  editorSelector: INPUT_EDITOR,
  parts: INPUT_PARTS,
});
kind('slider handle', {
  aliases: ['handle', 'thumb'],
  selector: '[role="slider"], svg[type="range-slider"] [name$="-handle"]',
  match: ['aria', 'name', 'text', 'dart'],
  dartNames: ['{q}-handle', '{q}'],
});

// --- forms ----------------------------------------------------------------------------------------
kind('property', {
  aliases: ['property row'],
  selector: '.u2-propgrid-row, tr.property-grid-item',
  match: ['label', 'name', 'aria', 'dart'],
  labelSelector: '.u2-propgrid-name, .property-grid-item-name-text',
  dartNames: ['prop-{q}'],
  parts: {label: '.u2-propgrid-name, .property-grid-item-name', value: '.u2-propgrid-value, .property-grid-item-value',
    editor: '[name^="prop-view-"]'},
});
// the Dart property grid's category is a row of the grid (prop-category-<name>)
kind('category', {
  aliases: ['property category'],
  selector: '.u2-propgrid-category, tr.property-grid-item[name^="prop-category-"]',
  match: ['title', 'dart', 'text'],
  labelSelector: '.u2-propgrid-category-title, .property-grid-item-name-text',
  dartNames: ['prop-category-{q}'],
});

// --- items of u2 collections: also addressed by key (data-u2-item), revealed for a gesture when not rendered
kind('item', {
  aliases: ['list item', 'row', 'option', 'entry'],
  selector: '.u2-list-row, [role="option"], [role="row"], [role="listitem"], .d4-list-item, li',
  match: ['text', 'label', 'aria', 'item', 'name'],
  items: true,
  labelSelector: `${PRIMARY_TEXT}, ${ITEM_LABEL}`,
});
// the Dart tree wraps every row in a .d4-tree-view-group[role=treeitem] that carries the row's
// text too, so the wrapper is excluded and the row — which owns name="tree-My-stuff---Favorites" —
// is the node
kind('tree node', {
  aliases: ['node', 'tree item'],
  selector: '[role="tree"] .u2-list-row, [role="treeitem"]:not(.d4-tree-view-group), .d4-tree-view-node',
  match: ['dart', 'label', 'text', 'item', 'name'],
  items: true,
  labelSelector: '.u2-tree-label, .d4-tree-view-node-label, .d4-tree-view-group-label, .d4-tree-view-item-label, ' + ITEM_LABEL,
  dartNames: ['tree-{q}'],
});
kind('table row', {
  selector: u2('table') + ' tbody tr, table tbody tr, [data-u2="data-table"] [data-u2-item]',
  match: ['label', 'text', 'item'],
  items: true,
  labelSelector: `td:first-child, th:first-child, ${ITEM_LABEL}`,
});

// --- display --------------------------------------------------------------------------------------
kind('heading', {aliases: ['header', 'title'], selector: 'h1, h2, h3, h4, h5, h6', match: ['text']});
kind('text', {aliases: ['label'], selector: '*', match: ['exact-text']});
kind('link', {selector: 'a, .d4-link-label', match: ['text', 'name', 'aria', 'dart'], dartNames: ['label-{q}']});
kind('action', {aliases: ['action link'], selector: '.d4-link-action', match: ['exact-text', 'text'],
  description: 'a function link of the context panel\'s Actions pane ("Convert Notation...")'});

// --- navigation -----------------------------------------------------------------------------------
// a group item (u2 nests its submenu inside the item, Dart its children) contains its
// children's labels, so an item's own label is the direct child
kind('menu item', {
  selector: '.u2-menu-item, [role="menuitem"], .d4-menu-item',
  match: ['label', 'text', 'item', 'name', 'dart'],
  items: true,
  labelSelector: `:scope > .u2-menu-label, :scope > .d4-menu-item-label, :scope > ${ITEM_LABEL}`,
  dartNames: ['div-{q}'],
});
// the platform's ribbon holds its commands as plain divs, named by the text or the icon they show
kind('ribbon item', {
  aliases: ['ribbon command'],
  selector: '.d4-ribbon-item',
  match: ['text', 'name', 'aria'],
});
kind('breadcrumb', {
  aliases: ['crumb'],
  selector: '.u2-breadcrumbs-item, .u2-breadcrumbs-current, ' + u2('breadcrumbs') + ' a',
  match: ['text', 'aria', 'item'],
  items: true,
});

// --- containers -----------------------------------------------------------------------------------
// .tab-handle is the Dart dock manager's tab, .dockspan-tab-handle the dock-spawn-ts one a function
// view keeps in its shadow root (which a CSS locator reaches); the label is a child either way
kind('tab', {
  selector: '[role="tab"], .d4-tab-header, .tab-handle, .dockspan-tab-handle',
  match: ['name', 'label', 'text', 'aria', 'item', 'dart'],
  items: true,
  labelSelector: '.tab-handle-text, .dockspan-tab-handle-text, ' + ITEM_LABEL,
  dartNames: ['{q}', 'tab-{q}'],
});
/** The tab of an open view in the shell's tab strip (`name="view-handle: <view>"`), by the view's
 * name; simple mode hides the strip, the status bar's Tabs toggle shows it. */
kind('view tab', {
  selector: '.tab-handle[name^="view-handle: "]',
  match: ['label', 'text'],
  labelSelector: '.tab-handle-text',
});
kind('tab panel', {aliases: ['tab page'], selector: '[role="tabpanel"], .d4-tab-content', match: ['name', 'aria']});
kind('accordion header', {
  aliases: ['pane header'],
  selector: u2('accordion') + ' [role="button"], .d4-accordion-pane-header',
  // a Dart pane header carries its own name; its text can hold a count too ("Activity2"), which no
  // exact-text match would find
  match: ['dart', 'text', 'aria'],
  dartNames: ['div-section--{q}'],
});
/** The entity cards of a Dart gallery (a connection's queries, the scripts browser, the Plugins
 * view's package cards): named `div-<Name>` and titled by their own label, not the
 * `.d4-item-card` the card kind knows. */
kind('gallery card', {
  selector: '.d4-gallery-card, .grok-app-card',
  match: ['dart', 'label', 'text'],
  // a package card's name block also holds its version
  labelSelector: '.grok-gallery-grid-item-title, .grok-package-name > div:first-child',
  dartNames: ['div-{q}'],
});
/** The toggles of the status bar's windows manager (Tabs, Toolbox, Context Panel, Console,
 * Presentation mode …), by their aria-label; `active` is the only state they publish. */
kind('status bar toggle', {selector: '.layout-status-bar [name^="toggle-"]', match: ['aria']});
kind('wizard step', {
  aliases: ['step'],
  selector: '.u2-wizard-step',
  match: ['label', 'text', 'item'],
  items: true,
  labelSelector: '.u2-wizard-title, ' + ITEM_LABEL,
});
kind('dock panel', {
  selector: '.panel-base',
  match: ['label', 'text'],
  labelSelector: '.panel-titlebar-text',
  description: 'a panel of the dock manager by the title its titlebar (or its tab handle) shows — "Activity cliffs"',
});
kind('splitter panel', {aliases: ['split pane'], selector: '.u2-splitter-panel', match: ['name']});
kind('sash', {aliases: ['splitter sash', 'divider'], selector: '.u2-sash, .u2-splitter-sash', match: ['aria', 'name']});

// --- the Dart shell -------------------------------------------------------------------------------
// the title bar belongs to the dock panel around the viewer, hence the ancestor walk
const PANEL = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]';
kind('viewer', {
  selector: '[name^="viewer-"], .d4-viewer',
  match: ['dart', 'name'],
  dartNames: ['viewer-{q}'],
  gestures: {click: 'mouse'},
  parts: {
    title: `${PANEL}//*[contains(@class, "panel-titlebar-text")]`,
    'settings icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="icon-font-icon-settings"]`,
    // the "?" of the title bar: its tooltip is the viewer's own summary above "Click for help (F1)"
    'help icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="icon-font-icon-help"]`,
    'menu icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="icon-font-icon-menu"]`,
    'close icon': `${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="Close" or @name="icon-font-icon-close"]`,
    description: '.d4-viewer-description',
    canvas: 'canvas[name="canvas"]',
    legend: '[name="legend"]',
    'mini legend icon': '[name="mini-legend-icon"]',
    'legend splitter': '[name="legend-splitter"]',
    'legend inner splitter': '[name="legend-inner-splitter"]',
    'legend close chevron': '[name="icon-hide-corner-legend"]',
    // the legend's one column selector, named after the property it is bound to, as every on-viewer selector
    'legend markers selector': '[name="legend"] [name^="div-column-combobox-"]',
    // the date map (year, quarter, …) the markers selector shows beside a datetime column
    'legend markers map': '[name="legend"] select[name^="input-aggr-selector-"]',
  },
});
// an entry of a viewer's legend: "R_ONE" legend item in legend of scatter plot viewer; selected while
// its category filters the viewer (aria-selected), named by the label it shows (aria-label)
kind('legend item', {
  aliases: ['legend entry'],
  selector: '[name="legend-item"]',
  match: ['aria', 'label', 'text'],
  labelSelector: '.d4-legend-value',
  parts: {label: '.d4-legend-value', cross: '.d4-legend-cross', thumbnail: 'canvas.d4-legend-value', marker: '[name="legend-item-marker"]'},
});
// a card of the filter panel by its caption: "RACE" filter card, checkbox of "RACE" filter card;
// aria-disabled on its body while suspended (the header, its checkbox and the mode word stay operable), its own counter as the indicator part
kind('filter card', {
  selector: '.d4-filter',
  match: ['label', 'dart'],
  labelSelector: '.d4-filter-column-name',
  dartNames: ['filter-card-{q}'],
  parts: {caption: '.d4-filter-column-name', checkbox: '.d4-filter-bool-input', indicator: '.d4-filter-indicator',
    mode: '[name="filter-mode-toggle"]', summary: '.d4-filter-summary', body: '.d4-filter-element',
    close: '[name="icon-times"]', 'search icon': '[name="icon-search"]', search: '.d4-filter-element input.d4-search-input'},
});
// the membership editor (membership_editor.dart): a member already in is a row, a search match a candidate
const MEMBERSHIP_PARTS = {checkbox: '.membership-row-admin input', 'checkbox label': '.membership-row-admin'};
kind('membership row', {
  aliases: ['member row'],
  selector: '.membership-row',
  match: ['label'],
  labelSelector: '.d4-user-selector-user-name',
  parts: {...MEMBERSHIP_PARTS, 'remove button': '[name="button-Remove"]'},
});
kind('membership candidate', {
  aliases: ['member candidate'],
  selector: '.membership-add-row',
  match: ['label'],
  labelSelector: '.d4-user-selector-user-name',
  parts: {...MEMBERSHIP_PARTS, 'add button': '[name="button-Add"]'},
});
// a table row of the Save project dialog's entity list (project_entity_move.dart): its action,
// its Data sync switch and the Creation script block
kind('project table', {
  aliases: ['saved table'],
  selector: '.grok-project-move-entity-row',
  match: ['dart'],
  dartNames: ['project-table-{q}'],
});
kind('view', {
  selector: '.d4-view-handle, [name^="view-handle: "]',
  match: ['dart', 'text'],
  dartNames: ['view-handle: {q}'],
});
kind('element', {
  aliases: ['component', 'control', 'panel', 'container', 'area', 'widget'],
  selector: '[data-u2-name], [data-widget]',
  match: ['name', 'dart'],
  dartNames: ['{q}'],
});
// a "name = value" status line: <div role="status"><label>name = </label><output>value</output></div>
kind('readout', {aliases: ['read-out'], selector: '[role="status"]', match: ['label', 'aria', 'name'],
  labelSelector: ':scope > label', parts: {label: ':scope > label', value: ':scope > output'},
  description: 'a "name = value" status line; "value of <name> readout" is the value alone'});
