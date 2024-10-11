/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {distinctUntilChanged, filter, take} from 'rxjs/operators';

import {ViewerApp as ViewerAppInstance} from './apps/ViewerApp';
import {FormApp as FormAppInstance} from './apps/FormApp';
import {HistoryApp as HistoryAppInstance} from './apps/HistoryApp';
import {ElementsApp as ElementsAppInstance} from './apps/ElementsApp';
import {TreeWizardApp as TreeWizardAppInstance} from './apps/TreeWizardApp';
import {SimpleDriverApp as SimpleDriverAppInstance} from './apps/SimpleDriverApp';
import {RFVWrapper} from './components/RFV/RFVWrapper';
import { PipelineConfiguration } from '@datagrok-libraries/compute-utils';

export const _package = new DG.Package();

//name: Tree Wizard
//tags: test, vue, model
//sidebar: @compute
//output: view app
//meta.icon: icons/tree-wizard.png
export async function TreeWizardApp() {
  return DG.Func.byName('Compute2:TreeWizardEditor').prepare({providerFunc: 'Compute2:MockProvider3'}).call();
}

//name: Tree Wizard Editor
//tags: editor
//input: string providerFunc
//output: view treeWizardView
export async function TreeWizardEditor(providerFunc: string) {
  const thisCall = grok.functions.getCurrentCall();

  await customElements.whenDefined('dg-markdown');

  const view = new DG.ViewBase();
  const app = Vue.createApp(TreeWizardAppInstance, {providerFunc});
  view.root.classList.remove('ui-panel');
  view.root.classList.add('ui-box');
  app.mount(view.root);
  
  view.name = 'DriverApp';
  view.parentCall = thisCall;
  view.parentView = thisCall.parentCall.aux['view'];
  view.basePath = `/${thisCall.func.name}`;
  grok.shell.addView(view);
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
