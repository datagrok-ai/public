/* Generic kinds: "<qualifier> <kind>" resolves without a registry entry. Every `data-u2` value the
   u2 library stamps (grep `dataset.u2 =` under libraries/u2/src) has a kind here — u2 is the
   standard UI, so the base vocabulary covers it whole; the Dart client's `name=` conventions
   (`annotate()` in html_utils.dart) ride along as `dartNames`, so the same phrase works on both UI
   generations. Inputs share the u2 parts (label, editor, options, error) as `<part> of <input>`. */
import {kind, KindDef, MatchStrategy} from '../../src/registry.js';

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

function u2(...ids: string[]): string {
  return ids.map((id) => `[data-u2="${id}"]`).join(', ');
}

function inputKind(name: string, ids: string[], dart: string, aliases: string[] = [], extra: Partial<KindDef> = {}): void {
  kind(name, {
    aliases,
    selector: [u2(...ids), dart].filter((s) => s).join(', '),
    match: INPUT_MATCH,
    labelSelector: INPUT_LABEL,
    dartNames: ['input-host-{q}'],
    editorSelector: INPUT_EDITOR,
    parts: INPUT_PARTS,
    ...extra,
  });
}

// --- inputs ---------------------------------------------------------------------------------------
kind('input', {
  aliases: ['field'],
  selector: '[data-u2$="-input"], ' + u2('combobox', 'multi-select', 'column-picker', 'column-combo', 'range-slider', 'typeahead') + ', .ui-input-root',
  match: INPUT_MATCH,
  labelSelector: INPUT_LABEL,
  dartNames: ['input-host-{q}'],
  editorSelector: INPUT_EDITOR,
  parts: INPUT_PARTS,
});
inputKind('text input', ['text-input'], '.ui-input-text', ['text field', 'textbox']);
inputKind('text area', ['text-area'], '.ui-input-textarea', ['textarea', 'multiline input']);
inputKind('choice input', ['choice-input'], '.ui-input-choice', ['dropdown', 'choice', 'select']);
inputKind('multi choice input', ['multi-choice-input'], '', ['multi choice']);
inputKind('number input', ['number-input', 'bigint-input', 'qnum-input'], '.ui-input-int, .ui-input-float', ['numeric input', 'number field']);
inputKind('checkbox', ['bool-input'], '.ui-input-bool, .u2-multi-choice-item, .u2-columns-option', ['bool input', 'switch', 'toggle'],
  {match: [...INPUT_MATCH, 'text']});
inputKind('date input', ['date-input', 'datetime-input'], '.ui-input-date', ['date field', 'datetime input', 'date picker']);
inputKind('color input', ['color-input'], '.ui-input-color', ['color picker']);
inputKind('font input', ['font-input'], '', ['font picker']);
inputKind('icon input', ['icon-input'], '', ['icon picker']);
inputKind('image input', ['image-input'], '', ['image picker']);
inputKind('list input', ['list-input'], '', []);
inputKind('map input', ['map-input'], '', ['key value input']);
inputKind('message input', ['message-input'], '', ['prompt input', 'chat input']);
inputKind('radio input', ['radio-input'], '.ui-input-radio', ['radio group', 'radio']);
inputKind('slider', ['slider-input'], '.ui-input-slider', ['slider input']);
inputKind('range slider', ['range-slider'], '', ['range input']);
kind('slider handle', {aliases: ['handle', 'thumb'], selector: '[role="slider"]', match: ['aria', 'name', 'text']});
inputKind('suggest input', ['suggest-input', 'typeahead'], '', ['autocomplete', 'typeahead']);
inputKind('combobox', ['combobox'], '.ui-input-editable-choice', ['combo box', 'editable choice']);
inputKind('multi select', ['multi-select'], '.ui-input-multi-choice', ['multi select input']);
inputKind('tags input', ['tags-input'], '.ui-input-tags', ['tags field']);
inputKind('file input', ['file-input', 'files-input'], '.ui-input-file', ['file field', 'files input']);
inputKind('columns input', ['columns-input', 'columns-map-input', 'aggregated-columns-input'], '.ui-input-columns', ['columns picker']);
inputKind('column input', ['column-combo', 'column-picker'], '.d4-column-selector', ['column picker', 'column selector', 'column combobox'],
  {match: ['name', 'label', 'dart'], dartNames: ['div-column-combobox-{q}-', 'input-host-{q}'], gestures: {open: 'mousedown'}});
inputKind('function input', ['function-input', 'func-call-input'], '', ['func input']);
inputKind('dynamic input', ['dynamic-input'], '', []);
inputKind('rsa input', ['rsa-input'], '', ['key input']);
inputKind('dart input', ['dart-input'], '', ['bridged input']);

// --- actions --------------------------------------------------------------------------------------
kind('button', {
  selector: 'button, [role="button"], .u2-btn, .ui-btn',
  match: ['name', 'text', 'aria', 'dart'],
  dartNames: ['button-{q}'],
});
kind('icon button', {
  selector: u2('icon-button') + ', button.u2-icon-btn',
  match: ['name', 'aria', 'dart'],
  dartNames: ['icon-{q}'],
});
kind('dropdown button', {
  aliases: ['menu button', 'split button'],
  selector: u2('dropdown-button'),
  match: ['name', 'text', 'aria'],
});
kind('button group', {aliases: ['segmented control'], selector: u2('button-group'), match: ['name', 'aria']});
kind('row actions', {aliases: ['actions'], selector: u2('row-actions'), match: ['name']});

// --- forms ----------------------------------------------------------------------------------------
kind('form', {
  selector: u2('form', 'object-form', 'func-form') + ', .ui-form',
  match: ['name', 'aria'],
});
kind('function form', {
  aliases: ['func form'],
  selector: u2('func-form'),
  match: ['name', 'aria'],
  parts: {'run button': u2('ff-run'), 'history icon': u2('ff-history-icon')},
});
kind('object form', {selector: u2('object-form'), match: ['name', 'aria']});
kind('property grid', {aliases: ['property editor'], selector: u2('property-grid', 'property-editor'), match: ['name', 'aria']});
kind('property', {aliases: ['property row'], selector: '.u2-propgrid-row', match: ['label', 'name'], labelSelector: '.u2-propgrid-name'});
kind('category', {
  aliases: ['property category'],
  selector: '.u2-propgrid-category',
  match: ['text', 'title'],
  labelSelector: '.u2-propgrid-category-title',
});

// --- collections ----------------------------------------------------------------------------------
kind('list', {
  selector: '[data-u2="list"]:not([role="tree"]), [role="listbox"], .u2-list:not([role="tree"]), .d4-list',
  match: ['name', 'aria'],
});
kind('item', {
  aliases: ['list item', 'row', 'option', 'entry'],
  selector: '.u2-list-row, [role="option"], [role="row"], .d4-list-item, li',
  match: ['text', 'label', 'aria', 'name'],
  labelSelector: PRIMARY_TEXT,
});
kind('tree', {selector: u2('tree') + ', [role="tree"], .d4-tree-view', match: ['name', 'aria']});
kind('tree node', {
  aliases: ['node', 'tree item'],
  selector: '[role="tree"] .u2-list-row, [role="treeitem"], .d4-tree-view-node',
  match: ['text', 'label', 'name'],
  labelSelector: '.u2-tree-label, .d4-tree-view-node-label',
});
kind('table', {selector: u2('table') + ', table', match: ['name', 'aria']});
kind('table row', {
  selector: u2('table') + ' tbody tr, table tbody tr',
  match: ['label', 'text'],
  labelSelector: 'td:first-child, th:first-child',
});
kind('grid', {aliases: ['virtual grid', 'cell grid'], selector: u2('grid') + ', .d4-grid', match: ['name', 'aria']});
kind('functions browser', {
  selector: u2('functions-browser'),
  match: ['name'],
  parts: {search: u2('fb-search'), list: u2('fb-list'), status: u2('fb-status'), 'clear button': u2('fb-clear'),
    roles: u2('fb-roles'), tags: u2('fb-tags'), 'empty message': u2('fb-empty')},
});
kind('history browser', {
  aliases: ['func call history browser', 'run history'],
  selector: u2('func-call-history-browser'),
  match: ['name'],
  parts: {list: u2('fch-list'), 'retry button': u2('fch-retry'), state: u2('fch-state')},
});

// --- display --------------------------------------------------------------------------------------
kind('icon', {
  selector: u2('icon') + ', .u2-icon, .grok-icon, [name^="icon-"], i[class*="fa-"]',
  match: ['aria', 'name', 'dart'],
  dartNames: ['icon-{q}'],
});
kind('badge', {selector: u2('badge', 'count-badge', 'dot'), match: ['text', 'name', 'aria']});
kind('tag', {selector: u2('tag') + ', .d4-tag', match: ['text', 'name']});
kind('notification', {
  aliases: ['balloon', 'toast', 'message'],
  selector: '[data-u2="notify"] .u2-notify, .d4-balloon',
  match: ['text'],
});
kind('progress bar', {
  aliases: ['progress'],
  selector: u2('progress-bar') + ', [role="progressbar"]',
  match: ['name', 'label', 'aria'],
  labelSelector: '.u2-progress-description',
});
kind('stat card', {aliases: ['stat'], selector: u2('stat-card'), match: ['name', 'label', 'text'], labelSelector: '.u2-stat-label'});
kind('tooltip', {selector: u2('tooltip') + ', [role="tooltip"], .d4-tooltip', match: ['text']});
kind('tour', {selector: u2('tour'), match: ['name']});
kind('async view', {aliases: ['loading view'], selector: u2('async-view'), match: ['name']});
kind('tray', {selector: u2('tray'), match: ['name']});
kind('heading', {aliases: ['header', 'title'], selector: 'h1, h2, h3, h4, h5, h6', match: ['text']});
kind('text', {aliases: ['label'], selector: '*', match: ['exact-text']});
kind('link', {selector: 'a, .d4-link-label', match: ['text', 'name', 'aria', 'dart'], dartNames: ['label-{q}']});

// --- navigation -----------------------------------------------------------------------------------
kind('toolbar', {selector: u2('toolbar') + ', [role="toolbar"], .d4-ribbon', match: ['name', 'aria']});
// the Dart main menu keeps a permanent [role="menu"] container per group — not a popup
kind('menu', {selector: u2('menu') + ', .d4-menu-popup, [role="menu"]:not(.d4-menu-item-container)', match: ['name', 'aria']});
kind('menu bar', {selector: u2('menu-bar') + ', [role="menubar"]', match: ['name', 'aria']});
kind('menu item', {
  selector: '.u2-menu-item, [role="menuitem"], .d4-menu-item',
  match: ['text', 'label', 'name', 'dart'],
  labelSelector: '.u2-menu-label, .d4-menu-item-label',
  dartNames: ['div-{q}'],
});
kind('breadcrumbs', {aliases: ['breadcrumb bar'], selector: u2('breadcrumbs'), match: ['name', 'aria']});
kind('breadcrumb', {
  aliases: ['crumb'],
  selector: '.u2-breadcrumbs-item, .u2-breadcrumbs-current, ' + u2('breadcrumbs') + ' a',
  match: ['text', 'aria'],
});

// --- containers -----------------------------------------------------------------------------------
kind('dialog', {
  selector: u2('dialog') + ', .d4-dialog, [role="dialog"]',
  match: ['name', 'title', 'dart'],
  labelSelector: '.u2-dialog-title-text, .d4-dialog-title',
  dartNames: ['dialog-{q}'],
  parts: {title: '.u2-dialog-title, .d4-dialog-title', 'close button': '.u2-dialog-close, .d4-dialog-close',
    footer: '.u2-dialog-footer, .d4-dialog-footer'},
});
kind('tabs', {aliases: ['tab strip', 'tab control'], selector: u2('tabs') + ', .d4-tab-control', match: ['name', 'aria']});
kind('tab', {
  selector: '[role="tab"], .d4-tab-header',
  match: ['name', 'text', 'aria', 'dart'],
  dartNames: ['{q}', 'tab-{q}'],
});
kind('tab panel', {aliases: ['tab page'], selector: '[role="tabpanel"], .d4-tab-content', match: ['name', 'aria']});
kind('section', {
  aliases: ['pane', 'accordion pane'],
  selector: u2('section') + ', .u2-section, .u2-accordion-pane, [name^="div-section--"], .d4-accordion-pane',
  match: ['name', 'title', 'dart'],
  labelSelector: '.u2-section-title, .u2-section-header, .u2-accordion-title, .d4-accordion-pane-header',
  dartNames: ['div-section--{q}', 'pane-{q}'],
});
kind('accordion', {selector: u2('accordion') + ', .d4-accordion', match: ['name', 'aria']});
kind('accordion header', {
  aliases: ['pane header'],
  selector: u2('accordion') + ' [role="button"], .d4-accordion-pane-header',
  match: ['text', 'aria'],
});
kind('card', {
  selector: u2('card', 'stat-card', 'entity-card') + ', .d4-item-card',
  match: ['name', 'title', 'text'],
  labelSelector: '.u2-card-title, .card-label',
});
kind('wizard', {selector: u2('wizard'), match: ['name', 'aria']});
kind('wizard step', {aliases: ['step'], selector: '.u2-wizard-step', match: ['label', 'text'], labelSelector: '.u2-wizard-title'});
kind('splitter', {aliases: ['split panel'], selector: u2('splitter'), match: ['name']});
kind('splitter panel', {aliases: ['split pane'], selector: '.u2-splitter-panel', match: ['name']});
kind('sash', {aliases: ['splitter sash', 'divider'], selector: '.u2-sash, .u2-splitter-sash', match: ['aria', 'name']});

// --- entities and the designer ---------------------------------------------------------------------
kind('chip', {selector: u2('chip', 'bind-chip', 'msg-chip'), match: ['text', 'name']});
kind('entity card', {aliases: ['entity'], selector: u2('entity-card'), match: ['name', 'text']});
kind('palette', {selector: u2('palette'), match: ['name']});
kind('designer', {selector: u2('designer'), match: ['name']});

// --- the Dart shell -------------------------------------------------------------------------------
kind('viewer', {
  selector: '[name^="viewer-"], .d4-viewer',
  match: ['dart', 'name'],
  dartNames: ['viewer-{q}'],
  gestures: {click: 'mouse'},
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
