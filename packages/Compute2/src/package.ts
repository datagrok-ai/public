/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {filter, take} from 'rxjs/operators';

import {ViewerApp as ViewerAppInstance} from './apps/ViewerApp';
import {FormApp as FormAppInstance} from './apps/FormApp';
import {HistoryApp as HistoryAppInstance} from './apps/HistoryApp';
import {ElementsApp as ElementsAppInstance} from './apps/ElementsApp';
import {TreeWizardApp as TreeWizardAppInstance} from './apps/TreeWizardApp';
import {SimpleDriverApp as SimpleDriverAppInstance} from './apps/SimpleDriverApp';
import {RFVApp} from './apps/RFVApp';
import { PipelineConfiguration } from '@datagrok-libraries/compute-utils';
import './tailwind.css'

export const _package = new DG.Package();

//tags: init
export async function init() {
  await DG.Func.byName('WebComponents:init').prepare().call();
}

//name: RichFunctionViewEditor
//tags: editor, vue
//input: funccall call
export async function RichFunctionViewEditor(call: DG.FuncCall) {
  const thisCall = grok.functions.getCurrentCall();

  await customElements.whenDefined('dg-markdown');

  const view = new DG.ViewBase();
  const app = Vue.createApp(RFVApp, {funcCall: call});
  grok.shell.add(view);
  app.mount(view.root);
  view.name = `${call.func.name}`;
  view.root.classList.remove('ui-panel');
  // view.root.classList.add('ui-box');
  view.root.style.overflow = 'hidden';

  view.name = call.func.friendlyName;
  view.parentCall = thisCall;
  view.parentView = thisCall.parentCall?.aux['view'];
  view.basePath = `/${call.func.name}`;

  grok.events.onViewRemoved.pipe(
    filter((closedView) => {
      return closedView === view;
    }),
    take(1),
  ).subscribe(() => {
    app.unmount();
  });
}


//name: Tree Wizard
//tags: test, vue, model
//sidebar: @compute
//meta.icon: icons/tree-wizard.png
export async function TreeWizardApp() {
  return DG.Func.byName('Compute2:TreeWizardEditor')
    .prepare({call: DG.Func.byName('Compute2:MockProvider3').prepare()}).call();
}

//name: Tree Wizard Editor
//tags: editor
//input: funccall call
export async function TreeWizardEditor(call: DG.FuncCall) {
  const thisCall = grok.functions.getCurrentCall();

  await customElements.whenDefined('dg-markdown');

  const view = new DG.ViewBase();
  const app = Vue.createApp(TreeWizardAppInstance, {providerFunc: call.func.nqName});
  view.root.classList.remove('ui-panel');
  view.root.classList.add('ui-box');
  grok.shell.add(view);
  app.mount(view.root);
  view.name = call.func.friendlyName;
  view.parentCall = thisCall;
  view.parentView = thisCall.parentCall?.aux['view'];
  view.basePath = `/${thisCall.func.name}`;

  grok.events.onViewRemoved.pipe(
    filter((closedView) => {
      return closedView === view;
    }),
    take(1),
  ).subscribe(() => {
    app.unmount();
  });
}

//tags: test, vue
export async function ViewerApp() {
  const view = new DG.ViewBase();
  const app = Vue.createApp(ViewerAppInstance);
  app.mount(view.root);
  view.name = 'ViewerApp';
  grok.shell.addView(view);
}

//tags: test, vue
export async function FormApp() {
  const view = new DG.ViewBase();
  const app = Vue.createApp(FormAppInstance);
  app.mount(view.root);
  view.name = 'FormApp';
  grok.shell.addView(view);
}


//tags: test, vue
export async function ElementsApp() {
  const view = new DG.ViewBase();
  const app = Vue.createApp(ElementsAppInstance);
  app.mount(view.root);
  view.name = 'ElementsApp';
  grok.shell.addView(view);
}

//tags: test, vue
export async function HistoryApp() {
  const view = new DG.ViewBase();
  const app = Vue.createApp(HistoryAppInstance);
  view.root.classList.remove('ui-panel');
  app.mount(view.root);
  view.name = 'HistoryApp';
  grok.shell.addView(view);
}

//tags: test, vue
export async function SimpleDriverApp() {
  // TODO: close view handling
  const view = new DG.ViewBase();
  const app = Vue.createApp(SimpleDriverAppInstance);
  view.root.classList.remove('ui-panel');
  app.mount(view.root);
  view.name = 'SimpleDriverApp';
  grok.shell.addView(view);
}


//name: MockWrapper1
export async function MockWrapper1() {}

//name: MockProvider1
//input: object params
//output: object result
export async function MockProvider1(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipeline1',
    nqName: 'Compute2:MockWrapper1',
    provider: 'Compute2:MockProvider1',
    version: '1.0',
    type: 'static',
    steps: [
      {
        id: 'step1',
        nqName: 'Compute2:LongScript',
      },
      // {
      //   id: 'step2',
      //   nqName: 'Compute2:LongFailingScript',
      // },
    ],
    links: [{
      id: 'link1',
      from: 'in1:step1/res',
      to: 'out1:step2/a',
    }],
  };
  return c;
}

//name: MockWrapper2
export async function MockWrapper2() {}


//name: MockProvider2
//input: object params
//output: object result
export async function MockProvider2(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipelinePar',
    nqName: 'Compute2:MockWrapper2',
    provider: 'Compute2:MockProvider2',
    version: '1.0',
    type: 'parallel',
    stepTypes: [{
      id: 'stepAdd',
      nqName: 'Compute:ObjectCooling',
      friendlyName: 'cooling',
    }, {
      id: 'stepMul',
      nqName: 'Compute2:TestMul2',
      friendlyName: 'mul',
    }, {
      type: 'ref',
      provider: 'Compute2:MockProvider1',
      version: '1.0',
    }],
    initialSteps: [
      {
        id: 'stepAdd',
      }, {
        id: 'pipeline1',
      },
    ],
  };
  return c;
}

//name: MockWrapper3
export async function MockWrapper3() {}

//name: MockProvider3
//input: object params
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockProvider3(params: any) {
  const c: PipelineConfiguration = {
    id: 'pipelinePar',
    nqName: 'Compute2:MockWrapper3', // for history
    provider: 'Compute2:MockProvider3', // for config
    friendlyName: 'Tree wizard model',
    version: '1.0',
    type: 'parallel',
    stepTypes: [{
      type: 'ref',
      provider: 'Compute2:MockProvider2',
      version: '1.0',
    }],
    initialSteps: [
      {
        id: 'pipelinePar',
      },
    ],
  };
  return c;
}

//name: TestAdd2
//input: double a
//input: double b
//output: double res
export async function TestAdd2(a: number, b: number) {
  return a + b;
}


//name: TestMul2
//input: double a
//input: double b
//output: double res
export async function TestSub2(a: number, b: number) {
  return a - b;
}

//input: double a
//input: double b
//output: double res
export async function TestMul2(a: number, b: number) {
  return a * b;
}

//input: double a
//input: double b
//output: double res
export async function TestDiv2(a: number, b: number) {
  return a / b;
}

//input: dataframe df
//output: dataframe res
export async function TestDF1(df: DG.DataFrame) {
  return df;
}
