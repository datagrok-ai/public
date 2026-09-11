/* The Datagrok shell by name. Loaded for every profile, so these names are reserved: an app
   registering "toolbox" for a bar of its own is refused — it registers "toolbar"-like names on its
   context instead, and "toolbox" keeps meaning the platform's. Selector sources: toolbox.dart
   (.d4-toolbox), console.dart, the shell's `name=` annotations (Browse, Toolbox), the selectors.ts
   files under playwright-public. */
import {element} from '../../src/registry.js';

element('toolbox', {selector: '.d4-toolbox', aliases: ['toolbox pane'],
  parts: {'viewers section': '[name="div-section--Viewers"]'}});
element('toolbox tab', {selector: '[name="Toolbox"]', aliases: ['toolbox sidebar tab']});
element('browse tab', {selector: '[name="Browse"]'});
element('browse panel', {selector: '.grok-view-browse, .layout-browse', aliases: ['browse view'],
  description: 'the left panel Browse opens; a bdd page starts in simple mode, where it is not in the DOM at all'});
element('browse tree', {selector: '.grok-view-browse [role="tree"], .layout-browse [role="tree"]',
  description: 'the tree inside the browse panel — the scope for a node phrase ("Files tree node inside browse tree")'});
element('context panel', {selector: '.grok-prop-panel', aliases: ['property panel']});
element('console', {selector: '.d4-console-wrapper'});
element('status bar', {selector: '.layout-status-bar', aliases: ['statusbar'],
  parts: {'view panel': '.d4-view-status-panel'}});
element('open tableview', {selector: '.d4-table-view, .grok-table-view', aliases: ['current table view', 'table view']});
element('grid', {selector: '[name="viewer-Grid"]', aliases: ['the grid'], gestures: {click: 'mouse'}});
element('gallery', {selector: '.grok-gallery-grid', aliases: ['item gallery'],
  description: 'the card gallery of the platform — the contents of a Files folder, a space, the Apps list'});
element('gallery search', {selector: '.grok-gallery-search-bar .ui-input-type-ahead'});
element('code editor', {selector: '.cm-editor, .CodeMirror', aliases: ['source editor'],
  description: 'the CodeMirror the platform embeds wherever code or a formula is edited — version 6 ' +
    'in the packages (.cm-editor), version 5 in the script view of the shell (.CodeMirror)'});
element('context menu', {selector: '.d4-menu-popup', aliases: ['popup menu'],
  description: 'the open Dart popup menu (the last one when a submenu is open)'});
element('cell editor', {selector: '[name="cell-editor"]', aliases: ['grid cell editor'],
  description: 'the value editor the grid shows over the cell being edited; absent between edits'});
element('filter panel', {selector: '[name="viewer-Filters"]', aliases: ['filters panel'],
  description: 'the Filters viewer of the current view (the same element as "filters viewer"); {widget} accepts it, so it has readings and hit areas of its own',
  parts: {counter: '[name="active-filter-counter"]', master: '[name="filters-master"]', search: '[name="filters-search"]',
    'add filter selector': '[name="div-column-combobox-add-filter"]', 'reset icon': '[name="icon-arrow-rotate-left"]',
    'search icon': '.d4-filter-group-header [name="icon-search"]', 'expand icon': '[name="icon-sort"]'}});
element('color picker icon', {selector: '[name="legend-icon-color-picker"]',
  description: 'the palette icon a hovered legend item shows to its left (the platform appends it to the page body)'});
element('marker picker icon', {selector: '[name="legend-icon-marker-picker"]',
  description: 'the shape icon a hovered marker item of a legend shows'});
