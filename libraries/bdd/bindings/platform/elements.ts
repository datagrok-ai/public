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
element('context panel', {selector: '.grok-prop-panel', aliases: ['property panel']});
element('console', {selector: '.d4-console-wrapper'});
element('status bar', {selector: '.layout-status-bar', aliases: ['statusbar'],
  parts: {'view panel': '.d4-view-status-panel'}});
element('open tableview', {selector: '.d4-table-view, .grok-table-view', aliases: ['current table view', 'table view']});
element('grid', {selector: '[name="viewer-Grid"]', aliases: ['the grid'], gestures: {click: 'mouse'}});
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
