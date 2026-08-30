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
import '@datagrok-libraries/u2/css/range-slider.css';
import '@datagrok-libraries/u2/css/radio.css';
import '@datagrok-libraries/u2/css/color.css';
import '@datagrok-libraries/u2/css/date.css';
import '@datagrok-libraries/u2/css/font.css';
import '@datagrok-libraries/u2/css/tags.css';
import '@datagrok-libraries/u2/css/choice.css';
import '@datagrok-libraries/u2/css/multi-select.css';
import '@datagrok-libraries/u2/css/combobox.css';
import '@datagrok-libraries/u2/css/list.css';
import '@datagrok-libraries/u2/css/functions-browser.css';
import '@datagrok-libraries/u2/css/grid.css';
import '@datagrok-libraries/u2/css/icon-input.css';
import '@datagrok-libraries/u2/css/function-input.css';
import '@datagrok-libraries/u2/css/func-call-history-browser.css';
import '@datagrok-libraries/u2/css/func-call-input.css';
import '@datagrok-libraries/u2/css/message-input.css';
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
import '@datagrok-libraries/u2/css/notify.css';
import '@datagrok-libraries/u2/css/tour.css';
import '@datagrok-libraries/u2/css/form.css';
import '@datagrok-libraries/u2/css/section.css';
import '@datagrok-libraries/u2/css/card.css';
import '@datagrok-libraries/u2/css/property-grid.css';
import '@datagrok-libraries/u2/css/async.css';
import '@datagrok-libraries/u2/css/typeahead.css';
import '@datagrok-libraries/u2/css/entity.css';
import '@datagrok-libraries/u2/css/badge.css';
import '@datagrok-libraries/u2/css/table.css';
import '@datagrok-libraries/u2/css/adaptive.css';
import '@datagrok-libraries/u2/css/buttons.css';
import '@datagrok-libraries/u2/css/designer.css';
import '@datagrok-libraries/u2/css/viewers.css';
import '../css/u2demo.css';

import {computed} from '@datagrok-libraries/u2';
import {appView, designerView, disposePanel, registerControlInspector, registerPlatformComponents,
  registerSpecNodeHandler} from '@datagrok-libraries/u2/src/dg/index.js';
import {DemoShell, buildDemo} from './demo';
import {DemoLeaf, leafForPath, pathOf} from './nav';
import {buildBrowseTree} from './browse';
import {demoRibbon} from './ribbon';
import {registerDemoSourceHandler} from './source-panel';
import {buildReportsBrowser} from './reports-browser';
import {registerEnabledEditors} from './editors';
import {registerPropRowHandler} from './convergence';
import {DESIGNER_SPEC, appendRunLog, designerContext} from './designer';

export * from './package.g';
// not a package function: what the e2e leak check reads after closing every view
export {leakReport, viewers} from '@datagrok-libraries/u2/src/dg/index.js';
export {sliceSymbol} from './source-panel';

export const _package = new DG.Package();

// JS apps own the whole path (the platform prefixes nothing); fragile: if U2Demo ever drops to a
// single app function, Func.url returns '' and the route collapses to '/apps/U2demo'
const APP_PATH = '/apps/U2demo/U2Demo';

let demoView: DG.ViewBase | null = null;
let demoShell: DemoShell | null = null;

// compared on `.dart`: the JS wrapper appView() built is not the object shell.views hands back
// for the same Dart view, so object identity never matches
grok.events.onViewRemoved.subscribe((v) => {
  if (demoView != null && v.dart === demoView.dart) {
    demoView = null;
    demoShell = null;
  }
});

function openDemo(leaf: DemoLeaf): void {
  // grok.shell.views is an Iterable, not an array (shell.ts:328)
  if (demoView != null && [...grok.shell.views].some((v) => v.dart === demoView!.dart)) {
    grok.shell.v = demoView;
    demoShell!.navigate(leaf.id);
    return;
  }
  grok.shell.addView(u2DemoApp(pathOf(leaf)));
}

//name: U2 Demo
//tags: app
//meta.browsePath: Dev
//input: string path { meta.url: true; optional: true }
//output: view result
export function u2DemoApp(path?: string): DG.ViewBase {
  registerEnabledEditors();
  registerPropRowHandler();
  registerControlInspector();
  registerDemoSourceHandler();
  const {content, shell, status} = buildDemo({initial: leafForPath(path)?.id});
  // the inspector's panel outlives this view otherwise, holding the control it last rendered
  content.own(disposePanel);
  const view = appView({
    name: 'U2 Demo', content, status,
    ribbon: [content.run(() => demoRibbon(shell))],
    path: computed(() => APP_PATH + pathOf(shell.current.value)),
  });
  demoView = view;
  demoShell = shell;
  return view;
}

//name: u2DemoTreeBrowser
//input: dynamic treeNode
//meta.role: appTreeBrowser
//meta.app: U2 Demo
export function u2DemoTreeBrowser(treeNode: DG.TreeViewGroup): void {
  buildBrowseTree(treeNode, openDemo);
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
  registerPlatformComponents();
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

/* What the designer's data sources are demonstrated against: a client-side table small enough to
   read whole, with a parameter worth binding to an input. */
const ORDERS = [
  {orderId: 1001, customer: 'Aspirin Labs', city: 'Kyiv', total: 1240, daysAgo: 2},
  {orderId: 1002, customer: 'Bayer', city: 'Lviv', total: 380, daysAgo: 5},
  {orderId: 1003, customer: 'Roche', city: 'Basel', total: 2150, daysAgo: 11},
  {orderId: 1004, customer: 'Novartis', city: 'Basel', total: 640, daysAgo: 24},
  {orderId: 1005, customer: 'Pfizer', city: 'New York', total: 1790, daysAgo: 45},
  {orderId: 1006, customer: 'Merck', city: 'Darmstadt', total: 920, daysAgo: 88},
];

//name: demoOrders
//description: Demo orders placed within the last N days — the data source demo of the u2 designer
//input: int days = 30
//output: dataframe orders
export function demoOrders(days: number): DG.DataFrame {
  const rows = ORDERS.filter((order) => order.daysAgo <= days);
  return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.create(0);
}

//name: u2Record
//description: Records a line in the U2 Designer's Run log — the function to wire a button to
//input: string text
//output: string entry
export function u2Record(text: string): string {
  const entry = `u2Record: ${text}`;
  grok.shell.info(entry);
  // the balloon is gone in five seconds — the Run log of the demo form is what stays
  appendRunLog(entry);
  return entry;
}
