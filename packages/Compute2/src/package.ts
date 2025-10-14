/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {filter, take} from 'rxjs/operators';

import {ViewerTestApp as ViewerAppInstance} from './apps/ViewerTestApp';
import {FormTestApp as FormAppInstance} from './apps/FormTestApp';
import {HistoryTestApp as HistoryAppInstance} from './apps/HistoryTestApp';
import {TreeWizardApp as TreeWizardAppInstance} from './apps/TreeWizardApp';
import {RFVApp} from './apps/RFVApp';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import './tailwind.css';
import {CustomFunctionView} from '@datagrok-libraries/compute-utils/function-views/src/custom-function-view';
import {HistoryApp} from './apps/HistoryApp';
import {Subject} from 'rxjs';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {deserialize, serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {OptimizerParams, runOptimizer} from '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';

declare global {
  var initialURLHandled: boolean;
}
declare let ENABLE_VUE_DEV_TOOLS: any;

export * from './package.g';
export const _package = new DG.Package();

function setViewHierarchyData(call: DG.FuncCall, view: DG.ViewBase) {
  view.parentCall = call.parentCall;

  if (view.parentCall?.aux?.view)
    view.parentView = view.parentCall.aux.view;

  if (call?.func?.name)
    view.basePath = `/${call.func.name}`;
}

function setVueAppOptions(app: Vue.App<any>) {
  app.config.compilerOptions.isCustomElement = (tag) => tag.startsWith('dg-') || tag === 'dock-spawn-ts';
  if (ENABLE_VUE_DEV_TOOLS)
    app.config.performance = true;
}

export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    try {
      await DG.Func.byName('WebComponents:init').prepare().call();
    } catch (e) {
      console.log(e);
      grok.shell.error(`WebComponents package init error`);
    }
  }

  @grok.decorators.editor({name: 'Custom Function View Editor'})
  static async CustomFunctionViewEditor(call: DG.FuncCall) : Promise<DG.ViewBase> {
    const view = (await call.call()).getOutputParamValue() as CustomFunctionView;
    setViewHierarchyData(call, view);

    await view.isReady.pipe(filter((x) => x), take(1)).toPromise();
    const updateFCBus = new Subject<DG.FuncCall>();

    const app = Vue.createApp(HistoryApp, {name: view.funcNqName, showHistory: view.showHistory, updateFCBus: Vue.markRaw(updateFCBus)});
    setVueAppOptions(app);

    const sub = updateFCBus.subscribe((fc) => {
      view.linkFunccall(fc);
      view.onAfterLoadRun(fc);
    });

    app.mount(view.historyRoot);

    grok.events.onViewRemoved.pipe(
      filter((closedView) => {
        return closedView === (view as any);
      }),
      take(1),
    ).subscribe(() => {
      sub.unsubscribe();
      app.unmount();
    });

    return view;
  }


  @grok.decorators.editor({name: 'Rich Function View Editor'})
  static async RichFunctionViewEditor(call: DG.FuncCall) : Promise<DG.ViewBase> {
    const view = new DG.ViewBase();
    setViewHierarchyData(call, view);

    const app = Vue.createApp(RFVApp, {funcCall: Vue.markRaw(call), view: Vue.markRaw(view)});
    view.root.classList.remove('ui-panel');
    view.root.classList.add('ui-box');
    setVueAppOptions(app);

    app.mount(view.root);

    grok.events.onViewRemoved.pipe(
      filter((closedView) => {
        return closedView === view;
      }),
      take(1),
    ).subscribe(() => {
      app.unmount();
    });

    grok.shell.windows.showHelp = false;

    return view;
  }


  @grok.decorators.editor({name: 'Tree Wizard Editor'})
  static async TreeWizardEditor(call: DG.FuncCall) : Promise<DG.ViewBase> {
    const providerFunc = call?.func?.options?.provider ?? call?.func?.nqName;
    if (!providerFunc)
      throw new Error(`Model ${call?.func?.name} has no provider`);

    const view = new DG.ViewBase();
    setViewHierarchyData(call, view);

    const modelName = call.options?.['title'] ?? call.func?.friendlyName ?? call.func?.name;
    const version = call.inputs.params?.version;
    let instanceConfig = deserialize(call.options?.instanceConfig ?? 'null');

    if (instanceConfig)
      instanceConfig = Vue.markRaw(instanceConfig);

    const {resolve} = call.aux;

    const app = Vue.createApp(TreeWizardAppInstance, {providerFunc, modelName, version, instanceConfig, resolve, view: Vue.markRaw(view)});
    view.root.classList.remove('ui-panel');
    view.root.classList.add('ui-box');
    setVueAppOptions(app);

    app.mount(view.root);

    grok.events.onViewRemoved.pipe(
      filter((closedView) => {
        return closedView === view;
      }),
      take(1),
    ).subscribe(() => {
      app.unmount();
      if (resolve)
        resolve();
    });

    grok.shell.windows.showHelp = false;

    return view;
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static async StartWorkflow(
    nqName: string,
    version: string,
    @grok.decorators.param({'type': 'object'}) instanceConfig?: PipelineInstanceConfig) {
    const func = DG.Func.byName(nqName);
    // @ts-ignore-next-line
    const {promise, resolve} = Promise.withResolvers();
    const call = func.prepare({version});
    call.options.instanceConfig = serialize(instanceConfig);
    call.aux.resolve = resolve;
    call.edit();
    return promise;
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static async RunOptimizer(
    @grok.decorators.param({'type': 'object'}) params: OptimizerParams,
  ) {
    const [,calls] = await runOptimizer(params);
    return calls;
  }

  @grok.decorators.func({
    tags: [
      'test',
      'vue',
    ],
  })
  static async ViewerTestApp() {
    const view = new DG.ViewBase();
    const app = Vue.createApp(ViewerAppInstance);
    setVueAppOptions(app);
    app.mount(view.root);
    view.name = 'ViewerTestApp';
    grok.shell.addView(view);
  }


  @grok.decorators.func({
    tags: [
      'test',
      'vue',
    ],
  })
  static async FormTestApp() {
    const view = new DG.ViewBase();
    const app = Vue.createApp(FormAppInstance);
    setVueAppOptions(app);
    app.mount(view.root);
    view.name = 'FormTestApp';
    grok.shell.addView(view);
  }


  @grok.decorators.func({
    tags: [
      'test',
      'vue',
    ],
  })
  static async HistoryTestApp() {
    const view = new DG.ViewBase();
    const app = Vue.createApp(HistoryAppInstance);
    view.root.classList.remove('ui-panel');
    setVueAppOptions(app);
    app.mount(view.root);
    view.name = 'HistoryTestApp';
    grok.shell.addView(view);
  }

  @grok.decorators.func({
    tags: [
      'test',
      'compute2',
    ],
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockPipeline1(
    @grok.decorators.param({type: 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'pipeline1',
      friendlyName: 'Pipeline 1',
      nqName: 'Compute2:MockPipeline1',
      version: '1.0',
      type: 'static',
      steps: [
        {
          id: 'step1',
          nqName: 'Compute2:LongScript',
        },
        {
          id: 'step2',
          nqName: 'Compute2:LongFailingScript',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }],
    };
    return c;
  }


  @grok.decorators.func({
    tags: [
      'test',
      'compute2',
    ],
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockPipeline2(
    @grok.decorators.param({'type': 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'pipelinePar',
      nqName: 'Compute2:MockPipeline2',
      friendlyName: 'Pipeline 2',
      version: '1.0',
      type: 'sequential',
      stepTypes: [{
        id: 'cooling',
        nqName: 'Compute2:ObjectCooling2',
        friendlyName: 'cooling',
        actions: [{
          id: 'action1',
          from: 'in:ambTemp',
          to: 'out:initTemp',
          position: 'none',
          handler({controller}) {
            controller.setAll('out', controller.getFirst('in') * 2);
            return;
          },
        },
        {
          id: 'action2',
          from: 'in:ambTemp',
          to: 'out:ambTemp',
          position: 'menu',
          menuCategory: 'Test',
          description: 'My menu action',
          handler({controller}) {
            controller.setAll('out', controller.getFirst('in') * 2);
            return;
          },
        },
        {
          id: 'action3',
          from: 'in:ambTemp',
          to: 'out:ambTemp',
          position: 'buttons',
          menuCategory: 'Test',
          description: 'My view action',
          handler({controller}) {
            controller.setAll('out', controller.getFirst('in') * 2);
            return;
          },
        }],
      }, {
        id: 'stepAdd',
        nqName: 'Compute2:TestAdd2',
        friendlyName: 'add',
      }, {
        id: 'stepMul',
        nqName: 'Compute2:TestMul2',
        friendlyName: 'mul',
      }, {
        id: 'LongScript',
        nqName: 'Compute2:LongScript',
        friendlyName: 'long',
      }, {
        type: 'ref',
        provider: 'Compute2:MockPipeline1',
        version: '1.0',
      }],
      initialSteps: [
        {
          id: 'stepAdd',
        }, {
          id: 'stepMul',
        }, {
          id: 'cooling',
        },
      ],
      links: [{
        id: 'selector',
        type: 'selector',
        from: 'in:cooling/ambTemp',
        to: ['out1:cooling/title', 'out2:cooling/description', 'out3:cooling/tags'],
        handler({controller}) {
          const val = controller.getFirst('in');
          controller.setDescriptionItem('out1', `Title ${val}`);
          controller.setDescriptionItem('out2', `Description ${val}`);
          controller.setDescriptionItem('out3', [`tag ${val}`]);
        },
      }, {
        id: 'link1',
        from: 'in1:stepAdd/res',
        to: 'out1:stepMul/a',
        type: 'meta',
        handler({controller}) {
          const addRes = controller.getFirst('in1');
          controller.setViewMeta('out1', {items: addRes > 0? ['0', '1', '2']: ['1', '2', '3']});
        },
      },
      {
        id: 'toMul',
        from: 'from:stepAdd/res',
        to: 'to:stepMul/b',
        defaultRestrictions: {to: 'restricted'},
      },
      {
        id: 'initialTempValidator',
        from: [
          'initTemp:cooling/initTemp',
          'ambTemp:cooling/ambTemp',
        ],
        to: [
          'toInitTemp:cooling/initTemp',
        ],
        actions: 'actions:cooling',
        type: 'validator',
        handler({controller}) {
          const initTemp = controller.getFirst('initTemp');
          const ambTemp = controller.getFirst('ambTemp');

          if (initTemp < ambTemp) {
            const action = controller.getValidationAction('actions', 'action1');

            if (!action) return;

            controller.setValidation('toInitTemp',
              ({errors: [
                {
                  description: `Initial temperature should be more than ambient temperature ${ambTemp}`,
                  actions: [{actionName: `Set reasonable initial temperature`, action}],
                },
              ]}));
          } else
            controller.setValidation('toInitTemp');
        },
      }],
    };
    return c;
  }


  @grok.decorators.func()
  static async TestAdd2(
    a: number,
    b: number): Promise<number> {
    return a + b;
  }


  @grok.decorators.func()
  static async TestSub2(
    a: number,
    b: number) : Promise<number> {
    return a - b;
  }


  @grok.decorators.func()
  static async TestMul2(
    a: number,
    b: number) : Promise<number> {
    return Number(a) * b;
  }


  @grok.decorators.func()
  static async TestDiv2(
    a: number,
    b: number) : Promise<number> {
    return a / b;
  }


  @grok.decorators.func()
  static async TestDF1(
    df: DG.DataFrame) : Promise<DG.DataFrame> {
    return df;
  }

  @grok.decorators.func({
    tags: [
      'test',
      'compute2',
    ],
    name: 'Custom View (Compute 2 Test)',
    editor: 'Compute2:CustomFunctionViewEditor',
  })
  static async TestCustomView() {
    const view = new MyView();
    return view;
  }
}


class MyView extends CustomFunctionView {
  private aIn?: DG.InputBase;
  private bIn?: DG.InputBase;
  private res?: DG.InputBase;

  constructor() {
    super('Compute2:TestAdd2');
    this.box = true;
  }

  override buildIO() {
    this.aIn = ui.input.float('a', {onValueChanged: (val) => this.funcCall!.inputs.a = val});
    this.bIn = ui.input.float('b', {onValueChanged: (val) => this.funcCall!.inputs.b = val});
    this.res = ui.input.float('res');
    this.res.enabled = false;
    return ui.div([
      this.aIn,
      this.bIn,
      ui.div([ui.bigButton('Run', async () => {
        await this.funcCall!.call();
        await this.saveRun(this.funcCall!);
        this.res!.value = this.funcCall?.outputs.res;
      })]),
      this.res,
    ]);
  }

  override async onAfterLoadRun(loadedRun: DG.FuncCall) {
    console.log(`Loaded run ${loadedRun.id}`);
    this.aIn!.value = this.funcCall!.inputs.a;
    this.bIn!.value = this.funcCall!.inputs.b;
    this.res!.value = this.funcCall?.outputs.res;
  }
}
