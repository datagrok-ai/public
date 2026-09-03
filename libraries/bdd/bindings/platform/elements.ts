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
