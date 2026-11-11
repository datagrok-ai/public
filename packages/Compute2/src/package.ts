/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BehaviorSubject, of} from 'rxjs';
import {concatMap, filter, finalize, map, take} from 'rxjs/operators';

import {ViewerTestApp as ViewerAppInstance} from './apps/ViewerTestApp';
import {FormTestApp as FormAppInstance} from './apps/FormTestApp';
import {TreeWizardApp as TreeWizardAppInstance} from './apps/TreeWizardApp';
import {RFVApp} from './apps/RFVApp';
import {provideDgViewService} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import type {PipelineConfiguration, ViewersHook} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import type {IRuntimePipelineMutationController} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/RuntimeControllers';
import './tailwind.css';
import {CustomFunctionView} from '@datagrok-libraries/compute-utils/function-views/src/custom-function-view';
import {attachRunStorage} from './components/RunStorage/custom-function-view-runs';
import {browseRuns} from './components/RunStorage/run-browser-dialog';
import {setVueAppOptions} from './utils';
import {openRunComparison} from './components/RunComparison/comparison-views';
import {normalizePipelineInstanceConfig} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import type {
  PipelineInstanceConfig, PipelineInstanceConfigInput,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {deserialize, serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {CONFIG_PATH} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {CUSTOM_VIEW_OPTION, loadViewersHook} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {applyCustomExport, callExcelReport} from './components/TreeWizard/tree-export';
import {callExports} from './run-exports';
import type {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import type {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {loadSavedTree, readRun, runWorkflow, saveWorkflowRun} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import type {RunState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import type {ItemMetadata} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/view/ViewCommunication';
import {runOptimizerFinalized} from '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';
import type {OptimizerParams} from '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';
import {ModelCatalogView,
  startModelCatalog,
  makeModelTreeBrowser,
  renderRestPanel,
  setModelCatalogEventHandlers,
  setModelCatalogHandler} from '@datagrok-libraries/compute-utils/model-catalog';
import dayjs from 'dayjs';

declare global {
  // eslint-disable-next-line no-var
  var initialURLHandled: boolean;
}

export * from './package.g';
import {_package, isRoleOnlyModelFilter} from './package-instance';
export {_package} from './package-instance';

let initRunned = false;
let startUriLoaded = false;

function setViewHierarchyData(call: DG.FuncCall, view: DG.ViewBase) {
  view.parentCall = call.parentCall;

  if (view.parentCall?.aux?.view)
    view.parentView = view.parentCall.aux.view;

  if (call?.func?.name) {
    const pcFunc = view.parentCall?.func;
    let prefix = '';
    if (pcFunc?.options?.role === 'app' && pcFunc.package) {
      const pkgUrl = pcFunc.package.meta?.url ?? `/${pcFunc.package.name}`;
      prefix = `/apps${pkgUrl}`;
    }
    view.basePath = `${prefix}/${call.func.name}`;
  }
}

// Compute2-only extension of the Model Hub view: adds the run comparison tool
// to the ribbon menu without touching the shared model-catalog code (used by Compute1).
class Compute2ModelCatalogView extends ModelCatalogView {
  constructor(viewName: string, roleOnlyFilter = false) {
    super(viewName, roleOnlyFilter);
    this.ribbonMenu.group('Tools')
      .item('Browse Runs...', () => grok.functions.call('Compute2:BrowseRuns'))
      .item('Compare Runs...', () => grok.functions.call('Compute2:CompareRuns'))
      .endGroup();
  }
}

const modelCatalogOptions = {
  _package,
  ViewClass: Compute2ModelCatalogView,
  segment: 'Modelhub',
  viewName: 'Model Hub',
  funcName: 'modelCatalog',
  setStartUriLoaded: () => startUriLoaded = true,
  getStartUriLoaded: () => startUriLoaded,
  roleOnlyFilter: false,
};

const providerOf = (nqName: string) => DG.Func.byName(nqName)?.options.provider ?? nqName;

export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    if (initRunned)
      return;
    initRunned = true;

    modelCatalogOptions.roleOnlyFilter = isRoleOnlyModelFilter();

    setModelCatalogHandler();
    setModelCatalogEventHandlers(modelCatalogOptions);

    try {
      await DG.Func.byName('WebComponents:init').prepare().call();
    } catch (e) {
      console.log(e);
      grok.shell.error(`WebComponents package init error`);
    }

    try {
      const compute1Init = DG.Func.byName('Compute:init');
      if (compute1Init)
        await compute1Init.prepare().call();
    } catch (e) {
      console.log(e);
      grok.shell.error(`Compute1 package init error`);
    }
  }


  @grok.decorators.func({name: 'renderRestPanel'})
  static async renderPanel(@grok.decorators.param({type: 'func'}) func: DG.Func) : Promise<DG.Widget> {
    return renderRestPanel(func as any) as any;
  }


  @grok.decorators.app({
    browsePath: 'Compute',
    name: 'Model Hub',
    outputs: [{type: 'view', name: 'result'}],
  })
  static modelCatalog() {
    const view = startModelCatalog(modelCatalogOptions);
    if (view && Array.from(grok.shell.views).includes(view)) {
      grok.shell.v = view;
      return null;
    }
    return view;
  }


  @grok.decorators.func({
    meta: {role: ' ', app: ' '},
  })
  static modelCatalogTreeBrowser(treeNode: DG.TreeViewGroup, browseView: DG.ViewBase) {
    makeModelTreeBrowser(treeNode as any, modelCatalogOptions.roleOnlyFilter);
  }


  @grok.decorators.editor({name: 'Custom Function View Editor'})
  static async CustomFunctionViewEditor(call: DG.FuncCall) : Promise<DG.ViewBase> {
    const view = (await call.call()).getOutputParamValue() as CustomFunctionView;
    setViewHierarchyData(call, view);

    await view.isReady.pipe(filter((x) => x), take(1)).toPromise();
    const detach = attachRunStorage(view, call.func);
    view.hostExports = callExports(call.func, true).map(({id, name}) => ({
      name, export: (fc: DG.FuncCall) => applyCustomExport(fc, id, {}, call.func),
    }));
    if (call.aux.initialRunId) {
      await view.loadRun(call.aux.initialRunId)
        .catch((e) => grok.shell.error(`Couldn't open the run: ${e?.message ?? e}`));
    }

    grok.events.onViewRemoved.pipe(
      filter((closedView) => {
        return closedView === (view as any);
      }),
      take(1),
    ).subscribe(() => detach());

    return view;
  }


  @grok.decorators.editor({name: 'Rich Function View Editor'})
  static async RichFunctionViewEditor(call: DG.FuncCall) : Promise<DG.ViewBase> {
    const view = DG.toJs(DG.toDart(new DG.ViewBase())) as DG.View;
    setViewHierarchyData(call, view);

    const viewersHook = await loadViewersHook(call.func);

    const app = Vue.createApp(RFVApp,
      {funcCall: Vue.markRaw(call), initialRunId: call.aux.initialRunId, viewersHook});
    view.root.classList.remove('ui-panel');
    view.root.classList.remove('ui-box');
    setVueAppOptions(app);
    provideDgViewService(app, view);

    app.mount(view.root);

    grok.events.onViewRemoved.pipe(
      filter((closedView) => {
        return closedView === view;
      }),
      take(1),
    ).subscribe(() => {
      app.unmount();
    });

    return view;
  }


  @grok.decorators.editor({name: 'Tree Wizard Editor'})
  static async TreeWizardEditor(call: DG.FuncCall) : Promise<DG.ViewBase> {
    const providerFunc = call?.func?.options?.provider ?? call?.func?.nqName;
    if (!providerFunc)
      throw new Error(`Model ${call?.func?.name} has no provider`);

    const view = DG.toJs(DG.toDart(new DG.ViewBase())) as DG.View;
    setViewHierarchyData(call, view);

    const modelName = call.options?.['title'] ?? call.func?.friendlyName ?? call.func?.name;
    const version = call.inputs.params?.version;
    let instanceConfig = deserialize(call.options?.instanceConfig ?? 'null');

    if (instanceConfig)
      instanceConfig = Vue.markRaw(instanceConfig);

    const {resolve, initialRunId} = call.aux;

    const app = Vue.createApp(TreeWizardAppInstance,
      {providerFunc, modelName, version, instanceConfig, resolve, initialRunId});
    view.root.classList.remove('ui-panel');
    view.root.classList.remove('ui-box');
    setVueAppOptions(app);
    provideDgViewService(app, view);

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

    return view;
  }


  @grok.decorators.func({
    name: 'Start Workflow',
    description: 'Launch a compute workflow (pipeline) by its qualified name and open its editor.',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async StartWorkflow(
    @grok.decorators.param({
      options: {description: 'Qualified name of the workflow function to launch'},
    }) nqName: string,
    @grok.decorators.param({options: {description: 'Workflow version to run'}}) version: string,
    @grok.decorators.param({
      'type': 'object',
      'options': {description: 'Optional initial pipeline configuration'},
    }) instanceConfig?: PipelineInstanceConfig) {
    const func = DG.Func.byName(nqName);
    // @ts-ignore-next-line
    const {promise, resolve} = Promise.withResolvers();
    const call = func.prepare({version});
    call.options.instanceConfig = serialize(instanceConfig);
    call.aux.resolve = resolve;
    call.edit();
    return promise;
  }


  @grok.decorators.func({
    name: 'RunWorkflow',
    description: 'Run a compute workflow without UI and return the state of every step.',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async RunWorkflow(
    @grok.decorators.param({options: {description: 'Qualified name of the workflow function to run'}})
      nqName: string,
    @grok.decorators.param({options: {description: 'Workflow version to run'}})
      version: string,
    @grok.decorators.param({type: 'object', options: {description: 'Optional initial pipeline configuration'}})
      instanceConfig?: PipelineInstanceConfigInput) {
    return runWorkflow(providerOf(nqName), version, instanceConfig && normalizePipelineInstanceConfig(instanceConfig));
  }


  @grok.decorators.func({
    name: 'SaveWorkflowRun',
    description: 'Save the state returned by RunWorkflow as a new workflow run in history.',
    outputs: [{type: 'funccall', name: 'run'}],
  })
  static async SaveWorkflowRun(
    @grok.decorators.param({options: {description: 'Qualified name of the workflow function'}})
      nqName: string,
    @grok.decorators.param({options: {description: 'Workflow version'}})
      version: string,
    @grok.decorators.param({type: 'object', options: {description: 'The state of a RunWorkflow result'}})
      state: RunState,
    @grok.decorators.param({type: 'object', options: {description: 'Optional title, description and tags'}})
      metadata?: ItemMetadata) {
    return saveWorkflowRun(providerOf(nqName), version, state, metadata);
  }


  @grok.decorators.func({
    name: 'ReportFuncCallExcel',
    description: 'The default Excel report of a run, as its export menu builds it: [Blob, ExcelJS.Workbook].',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async ReportFuncCallExcel(
    @grok.decorators.param({type: 'funccall', options: {description: 'The run to report'}}) funcCall: DG.FuncCall,
    @grok.decorators.param({type: 'object', options: {description: 'Validation results by input', optional: true}})
      validationState?: Record<string, ValidationResult>,
    @grok.decorators.param({type: 'object', options: {description: 'Consistency by input', optional: true}})
      consistencyState?: Record<string, ConsistencyInfo>,
    @grok.decorators.param({type: 'object', options: {
      description: 'The viewers hook of a workflow step, in place of the function\'s own', optional: true}})
      viewersHook?: ViewersHook,
    @grok.decorators.param({type: 'object', options: {description: 'Meta values by output, for the hook', optional: true}})
      metaState?: Record<string, BehaviorSubject<any>>) {
    return callExcelReport(funcCall, {
      viewersHook: viewersHook ?? await loadViewersHook(funcCall.func), metaState,
      validation: validationState, consistency: consistencyState,
    });
  }


  @grok.decorators.func({
    name: 'LoadWorkflowRun',
    description: 'Load a saved workflow run without UI and return the state of every step, as RunWorkflow does.',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async LoadWorkflowRun(
    @grok.decorators.param({options: {description: 'Meta FuncCall id of the saved workflow run'}}) id: string) {
    return loadSavedTree(id).pipe(
      concatMap(({tree}) => of(tree).pipe(map(readRun), finalize(() => tree.close()))),
    ).toPromise();
  }


  @grok.decorators.func({
    name: 'OpenWorkflowRun',
    // eslint-disable-next-line max-len
    description: 'Open a saved run by its FuncCall id — a workflow run in the Tree Wizard, or a single function run in its editor.',
  })
  static async OpenWorkflowRun(
    @grok.decorators.param({options: {description: 'Meta FuncCall id of the saved workflow run'}}) id: string) {
    const metaCall = await historyUtils.shallowLoadRun(id);
    // a run saved from a custom view opens in it again, unless that view is gone
    const customView = metaCall.options[CUSTOM_VIEW_OPTION];
    const call = ((customView && DG.Func.byName(customView)) || metaCall.func).prepare({});
    call.aux.initialRunId = id;
    if (`${call.func.options['editor'] ?? ''}`.startsWith('Compute2:')) {
      call.edit();
      return;
    }
    // workflow steps and nested pipelines are usually functions without an editor of their own
    const editor = metaCall.options[CONFIG_PATH] ? 'Compute2:TreeWizardEditor' : 'Compute2:RichFunctionViewEditor';
    const editorCall = await DG.Func.byName(editor).prepare({call}).call();
    grok.shell.addView(editorCall.getOutputParamValue());
  }


  @grok.decorators.func({
    name: 'Run Optimizer',
    description: 'Run parameter optimization (fitting) for a model and return the resulting function calls.',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async RunOptimizer(
    @grok.decorators.param({
      'type': 'object',
      'options': {description: 'Optimizer parameters: target function, variables, and objective'},
    }) params: OptimizerParams,
  ) {
    const fin = await runOptimizerFinalized(params);
    return fin.calls;
  }


  @grok.decorators.func({
    name: 'Compare Runs',
    description: 'Compare data across model runs: scalars or a single table column',
  })
  static async CompareRuns() {
    openRunComparison();
  }

  @grok.decorators.func({
    name: 'Browse Runs',
    description: 'Browse saved runs of any model in My runs and in spaces, and open or compare them',
  })
  static async BrowseRuns() {
    await browseRuns();
  }


  // Code for testing

  @grok.decorators.func()
  static async ViewerTestApp() {
    const view = new DG.ViewBase();
    const app = Vue.createApp(ViewerAppInstance);
    setVueAppOptions(app);
    app.mount(view.root);
    view.name = 'ViewerTestApp';
    grok.shell.addView(view);
  }


  @grok.decorators.func()
  static async FormTestApp() {
    const view = new DG.ViewBase();
    const app = Vue.createApp(FormAppInstance);
    setVueAppOptions(app);
    app.mount(view.root);
    view.name = 'FormTestApp';
    grok.shell.addView(view);
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'hook'}]})
  static async FailingViewersHookMaker() {
    throw new Error('broken hook');
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'hook'}]})
  static async SlowViewersHookMaker() {
    await new Promise((resolve) => setTimeout(resolve, 300));
    grok.events.fireCustomEvent('compute2-test-viewers-hook-made', null);
    return (ioName: string, type: string) =>
      grok.events.fireCustomEvent('compute2-test-viewers-hook-applied', {ioName, type});
  }


  @grok.decorators.func({
    name: 'Mock Pipeline 1',
    description: 'Sample static two-step workflow configuration used for testing the workflow engine.',
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
          nqName: 'Compute2:FastScript',
        },
      ],
      links: [{
        id: 'link1',
        from: 'in1:step1/res',
        to: 'out1:step2/a',
      }, {
        id: 'rootTitle',
        type: 'nodemeta',
        from: 'in:step1/a',
        to: 'title:title',
        handler({controller}) {
          controller.setDescriptionItem('title', 'Pipeline 1 title');
        },
      }],
    };
    return c;
  }


  @grok.decorators.func({
    name: 'Mock Pipeline 2',
    description: 'Sample sequential workflow configuration with links, actions, and validators for testing.',
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


  @grok.decorators.func({
    name: 'Mock Keyboard Pipeline',
    description: 'Workflow with fixed steps, a reorderable list and an invalid step, used for testing keyboard use.',
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockKeyboardPipeline(
    @grok.decorators.param({'type': 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'keyboardRoot',
      nqName: 'Compute2:MockKeyboardPipeline',
      friendlyName: 'Keyboard',
      version: '1.0',
      type: 'static',
      steps: [{
        id: 'fixedAdd',
        nqName: 'Compute2:TestAdd2',
        friendlyName: 'fixed add',
        initialValues: {a: -1, b: 1},
      }, {
        id: 'fixedMul',
        nqName: 'Compute2:TestMul2',
        friendlyName: 'fixed mul',
        initialValues: {a: 2, b: 3},
      }, {
        id: 'list',
        type: 'dynamic',
        friendlyName: 'list',
        stepTypes: [{
          id: 'item',
          nqName: 'Compute2:TestAdd2',
          friendlyName: 'item',
        }, {
          id: 'pinned',
          nqName: 'Compute2:TestMul2',
          friendlyName: 'pinned',
          disableUIDragging: true,
          disableUIRemoving: true,
        }],
        initialSteps: [
          {id: 'item', initialValues: {a: 1, b: 0}},
          {id: 'item', initialValues: {a: 2, b: 0}},
          {id: 'item', initialValues: {a: 3, b: 0}},
          {id: 'pinned', initialValues: {a: 1, b: 1}},
        ],
        links: [{
          id: 'itemTitle',
          type: 'nodemeta',
          base: 'base:expand(item)',
          from: 'a:same(@base, item)/a',
          to: 'title:same(@base, item)/title',
          handler({controller}) {
            controller.setDescriptionItem('title', `item ${controller.getFirst('a')}`);
          },
        }],
      }, {
        id: 'slow',
        nqName: 'Compute2:LongScript',
        friendlyName: 'slow',
        initialValues: {a: 1, b: 1},
      }, {
        id: 'group',
        type: 'static',
        friendlyName: 'group',
        steps: [{id: 'grouped', nqName: 'Compute2:TestAdd2', friendlyName: 'grouped step'}],
      }],
      links: [{
        id: 'firstItemB',
        from: 'in:fixedMul/a',
        to: 'out:list/first(item)/b',
        defaultRestrictions: {out: 'restricted'},
      }, {
        id: 'fixedAddPositive',
        type: 'validator',
        from: 'a:fixedAdd/a',
        to: 'out:fixedAdd/a',
        handler({controller}) {
          controller.setValidation('out',
            controller.getFirst('a') < 0 ? {errors: ['a must not be negative']} : undefined);
        },
      }],
    };
    return c;
  }


  @grok.decorators.func({
    name: 'Mock Tree Wizard Pipeline',
    description: 'Workflow with step actions, a nested workflow, an action step, step history, a custom export ' +
      'and a return, used for testing the Tree Wizard.',
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockTreeWizardPipeline(
    @grok.decorators.param({'type': 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'treeWizardRoot',
      nqName: 'Compute2:MockTreeWizardPipeline',
      friendlyName: 'Tree wizard',
      version: '1.0',
      type: 'static',
      steps: [{
        id: 'first',
        nqName: 'Compute2:TestAdd2',
        friendlyName: 'first',
        initialValues: {a: 3, b: 1},
        enableHistory: true,
        actions: [{
          id: 'double',
          friendlyName: 'Double a',
          position: 'buttons',
          from: 'in:a',
          to: 'out:a',
          handler({controller}) {
            controller.setAll('out', controller.getFirst('in') * 2);
          },
        }, {
          id: 'reset',
          friendlyName: 'Reset a',
          position: 'buttons',
          confirmationMessage: 'Reset a to 1?',
          from: [],
          to: 'out:a',
          handler({controller}) {
            controller.setAll('out', 1);
          },
        }],
      }, {
        type: 'ref',
        provider: 'Compute2:MockTreeWizardGroup',
        version: '1.0',
      }, {
        id: 'review',
        type: 'action',
        friendlyName: 'review',
      }, {
        id: 'last',
        nqName: 'Compute2:TestAdd2',
        friendlyName: 'last',
        initialValues: {a: 1, b: 1},
      }],
      links: [{
        id: 'rootTitle',
        type: 'nodemeta',
        from: 'in:first/a',
        to: 'title:title',
        handler({controller}) {
          controller.setDescriptionItem('title', 'Wizard title');
        },
      }, {
        id: 'firstTitle',
        type: 'nodemeta',
        from: 'a:first/a',
        to: 'title:first/title',
        handler({controller}) {
          controller.setDescriptionItem('title', `first ${controller.getFirst('a')}`);
        },
      }],
      onReturn: {
        id: 'returnA',
        type: 'return',
        from: 'a:first/a',
        to: [],
        handler({controller}) {
          controller.returnResult({a: controller.getFirst('a')});
        },
      },
      customExports: [{
        id: 'recorder',
        friendlyName: 'Record export',
        handler: async (state, utils) => DG.FileInfo.fromString('recorder.txt',
          `Exported ${state.configId} with ${Object.keys(utils).sort().join(', ')}`),
      }],
    };
    return c;
  }


  @grok.decorators.func({
    name: 'Mock Tree Wizard Group',
    description: 'Nested workflow of the Mock Tree Wizard Pipeline.',
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockTreeWizardGroup(
    @grok.decorators.param({'type': 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'treeWizardGroup',
      nqName: 'Compute2:MockTreeWizardGroup',
      friendlyName: 'group',
      version: '1.0',
      type: 'static',
      steps: [{
        id: 'inner',
        nqName: 'Compute2:TestMul2',
        friendlyName: 'inner',
        initialValues: {a: 2, b: 3},
      }],
      links: [{
        id: 'groupTitle',
        type: 'nodemeta',
        from: 'in:inner/a',
        to: 'title:title',
        handler({controller}) {
          controller.setDescriptionItem('title', 'group title');
        },
      }],
    };
    return c;
  }


  @grok.decorators.func({
    name: 'Mock Single Step Pipeline',
    description: 'Single-step workflow used for testing the compact Tree Wizard view.',
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockSingleStepPipeline(
    @grok.decorators.param({type: 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'singleStep',
      friendlyName: 'Single Step',
      nqName: 'Compute2:MockSingleStepPipeline',
      version: '1.0',
      type: 'static',
      compactView: true,
      steps: [{
        id: 'cooling',
        nqName: 'Compute2:ObjectCooling2',
        friendlyName: 'Cooling',
        initialValues: {ambTemp: 20},
        inputRestrictions: {ambTemp: 'restricted'},
        actions: [{
          id: 'doubleInitTemp',
          from: 'in:initTemp',
          to: 'out:initTemp',
          position: 'buttons',
          friendlyName: 'Double initial temperature',
          handler({controller}) {
            controller.setAll('out', controller.getFirst('in') * 2);
          },
        }, {
          id: 'resetSimTime',
          from: 'in:simTime',
          to: 'out:simTime',
          position: 'menu',
          menuCategory: 'Test',
          friendlyName: 'Reset simulation time',
          handler({controller}) {
            controller.setAll('out', 3600);
          },
        }],
      }],
      links: [{
        id: 'initTempValidator',
        type: 'validator',
        from: ['initTemp:cooling/initTemp', 'ambTemp:cooling/ambTemp'],
        to: 'toInitTemp:cooling/initTemp',
        handler({controller}) {
          const tooCold = controller.getFirst('initTemp') < controller.getFirst('ambTemp');
          controller.setValidation('toInitTemp',
            tooCold ? {errors: [{description: 'Initial temperature should be above ambient'}]} : undefined);
        },
      }, {
        id: 'hideArea',
        type: 'meta',
        from: 'desiredTemp:cooling/desiredTemp',
        to: 'area:cooling/area',
        handler({controller}) {
          controller.setViewMeta('area', {hidden: controller.getFirst('desiredTemp') > 1000});
        },
      }],
    };
    return c;
  }


  @grok.decorators.func({
    name: 'Mock Single Step Nested',
    description: 'Single-step workflow reached through a ref, used for testing the compact Tree Wizard view.',
    editor: 'Compute2:TreeWizardEditor',
    outputs: [{type: 'object', name: 'result'}],
  })
  static async MockSingleStepNested(
    @grok.decorators.param({type: 'object'}) params: any) {
    const c: PipelineConfiguration = {
      id: 'singleStepNested',
      friendlyName: 'Single Step Nested',
      nqName: 'Compute2:MockSingleStepNested',
      version: '1.0',
      type: 'static',
      compactView: true,
      steps: [{
        type: 'ref',
        provider: 'Compute2:MockSingleStepPipeline',
        version: '1.0',
      }],
    };
    return c;
  }


  @grok.decorators.func({
    editor: 'Compute2:TreeWizardEditor',
    tags: ['stress'],
    outputs: [{type: 'object', name: 'result'}],
  })
  static async StressTestPipeline(
    @grok.decorators.param({type: 'object'}) params: any) {
    const PIPELINE_COUNT = 10;
    const STEPS_PER_PIPELINE = 20;

    const dynPipeConfig = {
      id: 'dynPipe',
      type: 'dynamic' as const,
      friendlyName: 'Pipeline',
      stepTypes: [{
        id: 'mutateStep',
        nqName: 'Compute2:StressMutateStep',
        friendlyName: 'Mutate',
        disableUIControlls: true,
      }],
      initialSteps: Array.from({length: STEPS_PER_PIPELINE}, () => ({id: 'mutateStep'})),
      onInit: {
        id: 'initFirstDf',
        from: [],
        to: 'out1:first(mutateStep)/df',
        handler({controller}: any) {
          const n = 100;
          const df = DG.DataFrame.create(n);
          df.name = 'Initial Data';
          df.columns.addNewInt('id').init((i: number) => i);
          df.columns.addNewFloat('val1').init(() => Math.random());
          df.columns.addNewFloat('val2').init(() => Math.random());
          df.columns.addNewString('cat1').init((i: number) => `cat_${i % 10}`);
          df.columns.addNewString('cat2').init((i: number) => `grp_${i % 5}`);
          df.columns.addNewBool('active').init((i: number) => i % 2 === 0);
          df.columns.addNewInt('count').init(() => Math.floor(Math.random() * 1000));
          df.columns.addNewFloat('score').init(() => Math.random() * 100);
          df.columns.addNewDateTime('timestamp').init(() => dayjs(Date.now() - Math.random() * 1e9));
          df.columns.addNewFloat('noise').init(() => Math.random());
          controller.setAll('out1', df);
        },
      },
      links: [
        {
          id: 'chain',
          base: 'base:expand(mutateStep)',
          from: 'from:same(@base, mutateStep)/dfOut',
          to: 'to:after+(@base, mutateStep)/df',
        },
      ],
    };

    const c: PipelineConfiguration = {
      id: 'stressRoot',
      friendlyName: 'Stress Test Pipeline',
      nqName: 'Compute2:StressTestPipeline',
      version: '1.0',
      type: 'dynamic',
      stepTypes: [
        {
          id: 'sharedProducer',
          nqName: 'Compute2:StressSharedProducer',
          friendlyName: 'Shared Data Producer',
          disableUIControlls: true,
        },
        {
          ...dynPipeConfig,
          disableUIControlls: true,
        },
      ],
      initialSteps: [
        {id: 'sharedProducer'},
        ...Array.from({length: PIPELINE_COUNT}, () => ({id: 'dynPipe'})),
      ],
      actions: [
        {
          id: 'addPipeline',
          type: 'pipeline' as const,
          position: 'buttons' as const,
          friendlyName: 'Add Pipeline',
          icon: 'plus',
          from: [],
          to: 'out1',
          handler({controller}: {controller: IRuntimePipelineMutationController}) {
            controller.addStep('out1', 'dynPipe');
          },
        },
      ],
      links: [
        {
          id: 'sharedLink',
          from: 'in1:sharedProducer/sharedDf',
          to: 'out1:all(dynPipe)/all(mutateStep)/sharedDf',
        },
      ],
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
    name: 'Custom View (Compute 2 Test)',
    editor: 'Compute2:CustomFunctionViewEditor',
    // meta: { role: 'model' },
    meta: {customExports: '[{"name":"rec","function":"Compute2:TestCustomExportRecorder"}]'},
    outputs: [{type: 'object', name: 'result'}],
  })
  static async TestCustomView() {
    const view = new MyView();
    return view;
  }


  // Fixtures for the custom-export report-handler test (see test/custom-export.ts).
  @grok.decorators.func({
    // eslint-disable-next-line max-len
    meta: {customExports: '[{"name":"rec","function":"Compute2:TestCustomExportRecorder"},{"name":"none","function":"Compute2:TestCustomExportNoFile"},{"name":"hook","function":"Compute2:TestCustomExportHookArgs"}]'},
  })
  static async TestCustomExportModel(a: number): Promise<number> {
    return a;
  }

  @grok.decorators.func()
  static async TestCustomExportRecorder(funcCall: DG.FuncCall): Promise<DG.FileInfo> {
    return DG.FileInfo.fromString('recorded.txt', `${funcCall?.func?.nqName}|${funcCall?.inputs?.['a']}`);
  }

  @grok.decorators.func()
  static async TestCustomExportNoFile(funcCall: DG.FuncCall): Promise<string> {
    return funcCall.func.nqName;
  }

  @grok.decorators.func()
  static async TestCustomExportHookArgs(
    funcCall: DG.FuncCall,
    @grok.decorators.param({type: 'object', options: {optional: true}}) viewersHook?: ViewersHook,
    @grok.decorators.param({type: 'object', options: {optional: true}}) metaState?: Record<string, BehaviorSubject<any>>,
  ): Promise<DG.FileInfo> {
    return DG.FileInfo.fromString('hook.txt', `${typeof viewersHook}|${metaState?.['a']?.value}`);
  }

  // Fixture for the RFV tab tests (see test/rfv-view.ts): empty outputs get no tab, so `last` has one only
  // while `rows` > 0
  @grok.decorators.func({outputs: [
    {type: 'dataframe', name: 'first'},
    {type: 'dataframe', name: 'middle'},
    {type: 'dataframe', name: 'last'},
  ]})
  static async TestRfvTabs(
    @grok.decorators.param({type: 'int'}) rows: number,
  ): Promise<any> {
    return {
      first: DG.DataFrame.fromCsv('x\n1'), middle: DG.DataFrame.fromCsv('y\n1'), last: DG.DataFrame.create(rows),
    };
  }

  // Fixtures for the project-export data-prep tests (see test/project-export-data.ts);
  // the tests only prepare() these calls, the bodies never run.
  @grok.decorators.func({outputs: [{type: 'dataframe', name: 'res'}]})
  static async TestProjectExportSingleOut(
    @grok.decorators.param({type: 'int'}) x: number,
  ): Promise<DG.DataFrame> {
    return DG.DataFrame.create(x);
  }

  @grok.decorators.func({outputs: [
    {type: 'dataframe', name: 'res1'},
    {type: 'dataframe', name: 'res2'},
  ]})
  static async TestProjectExportMultiOut(
    df: DG.DataFrame,
    @grok.decorators.param({type: 'int'}) x: number,
  ): Promise<any> {
    return {res1: df, res2: df};
  }

  @grok.decorators.func({outputs: [{type: 'dataframe', name: 'res'}]})
  static async TestProjectExportFileIn(
    @grok.decorators.param({options: {optional: true, nullable: true}}) file: DG.FileInfo,
    @grok.decorators.param({type: 'int'}) x: number,
  ): Promise<DG.DataFrame> {
    return DG.DataFrame.create(x);
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
        await this.funcCall!.call(undefined, undefined, {processed: true, report: false});
        this.funcCall!.started = dayjs();
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
