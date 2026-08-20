/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/* Every skin the library ships except css/icons.css (the standalone `.u2-standalone` icon font —
   the platform already serves Font Awesome) and css/theme-dark.css (an opt-in theme). A control
   whose skin is missing falls back to UA chrome, which reads as a u2 defect. */
import '@datagrok-libraries/u2/css/tokens.css';
import '@datagrok-libraries/u2/css/elements.css';
import '@datagrok-libraries/u2/css/inputs.css';
import '@datagrok-libraries/u2/css/number.css';
import '@datagrok-libraries/u2/css/slider.css';
import '@datagrok-libraries/u2/css/radio.css';
import '@datagrok-libraries/u2/css/color.css';
import '@datagrok-libraries/u2/css/date.css';
import '@datagrok-libraries/u2/css/font.css';
import '@datagrok-libraries/u2/css/tags.css';
import '@datagrok-libraries/u2/css/choice.css';
import '@datagrok-libraries/u2/css/multi-select.css';
import '@datagrok-libraries/u2/css/combobox.css';
import '@datagrok-libraries/u2/css/list.css';
import '@datagrok-libraries/u2/css/file.css';
import '@datagrok-libraries/u2/css/progress.css';
import '@datagrok-libraries/u2/css/spec.css';
import '@datagrok-libraries/u2/css/wizard.css';
import '@datagrok-libraries/u2/css/tree.css';
import '@datagrok-libraries/u2/css/tabs.css';
import '@datagrok-libraries/u2/css/splitter.css';
import '@datagrok-libraries/u2/css/accordion.css';
import '@datagrok-libraries/u2/css/breadcrumbs.css';
import '@datagrok-libraries/u2/css/toolbar.css';
import '@datagrok-libraries/u2/css/menu.css';
import '@datagrok-libraries/u2/css/menu-bar.css';
import '@datagrok-libraries/u2/css/dialog.css';
import '@datagrok-libraries/u2/css/tooltip.css';
import '@datagrok-libraries/u2/css/form.css';
import '@datagrok-libraries/u2/css/property-grid.css';
import '@datagrok-libraries/u2/css/async.css';
import '@datagrok-libraries/u2/css/typeahead.css';
import '@datagrok-libraries/u2/css/entity.css';
import '@datagrok-libraries/u2/css/badge.css';
import '@datagrok-libraries/u2/css/table.css';
import '@datagrok-libraries/u2/css/adaptive.css';
import '@datagrok-libraries/u2/css/buttons.css';
import '@datagrok-libraries/u2/css/designer.css';
import '../css/u2demo.css';

import {MenuBar, registerAll, signal} from '@datagrok-libraries/u2';
import {appView, designerView, registerSpecNodeHandler} from '@datagrok-libraries/u2/src/dg/index.js';
import {buildDemo} from './demo';
import {buildReportsBrowser} from './reports-browser';
import {registerEnabledEditors} from './editors';
import {registerPropRowHandler} from './convergence';
import {DESIGNER_SPEC, designerContext} from './designer';

export * from './package.g';

export const _package = new DG.Package();

//name: U2 Demo
//tags: app
//meta.browsePath: Dev
//output: view result
export function u2DemoApp(): DG.ViewBase {
  registerEnabledEditors();
  registerPropRowHandler();
  const autoRefresh = signal(false);
  const bar = new MenuBar()
    .menu('File', (m) => m
      .item('New session', () => grok.shell.info('New session'), {shortcut: 'Ctrl+N'})
      .item('Save layout', () => grok.shell.info('Layout saved'))
      .separator()
      .group('Export')
      .item('As JSON', () => grok.shell.info('Exported as JSON'))
      .item('As CSV', () => grok.shell.info('Exported as CSV'))
      .endGroup())
    .menu('View', (m) => m
      .item('Auto-refresh', () => autoRefresh.value = !autoRefresh.value,
        {check: autoRefresh.value})
      .item('Reset panels', () => grok.shell.info('Panels reset')))
    .menu('Help', (m) => m
      .item('About u2', () => grok.shell.info('u2: next-gen Datagrok UI library'))
      .item('Sources', () => window.open('https://github.com/datagrok-ai/public/tree/master/libraries/u2')));
  return appView({name: 'U2 Demo', content: buildDemo(), ribbon: [[bar]]});
}

//name: Reports Browser
//tags: app
//meta.browsePath: Dev
//output: view result
export function reportsBrowserApp(): DG.ViewBase {
  const {content, ribbon, status} = buildReportsBrowser();
  return appView({name: 'Reports', content, ribbon, status});
}

//name: U2 Designer
//tags: app
//meta.browsePath: Dev
//output: view result
export function u2DesignerApp(): DG.ViewBase {
  registerAll();
  registerSpecNodeHandler();
  return designerView(DESIGNER_SPEC, {name: 'U2 Designer', ctx: designerContext()});
}

//name: u2AutoRegisterEditors
//tags: autostart
//description: Registers u2 value editors for the property types enabled in `u2.valueEditors`
export function u2AutoRegisterEditors(): void {
  registerEnabledEditors();
}

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}
