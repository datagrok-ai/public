import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {DBSchema} from 'idb';
import {
  Viewer, InputForm,
  BigButton, IconFA,
  RibbonPanel, DockManager, MarkDown,
  RibbonMenu,
  ifOverlapping,
  tooltip,
  Button,
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FittingView, TargetDescription} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {SensitivityAnalysisView} from '@datagrok-libraries/compute-utils';
import {RangeDescription} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {ScalarsPanel, ScalarState} from './ScalarsPanel';
import {BehaviorSubject} from 'rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {useViewersHook} from '../../composables/use-viewers-hook';
import {ViewAction} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {useLayoutDb} from '../../composables/use-layout-db';
import {take} from 'rxjs/operators';

type PanelsState = {
  historyHidden: boolean,
  helpHidden: boolean,
  formHidden: boolean,
  visibleTabLabels: string[],
  layout: string,
};

interface ScalarsState {
  type: 'scalars',
  scalarsData: ScalarState[],
}

interface DataFrameState {
  name: string,
  df: DG.DataFrame,
  type: 'dataframe',
  config: Record<string,any>,
}

type TabContent = Map<string, ScalarsState | DataFrameState>;

interface RenderStateItem {
  tabLabel: string;
  tabContent: ScalarsState | DataFrameState;
  isInput: boolean;
}

const getEmptyTabToProperties = () => ({
  inputs: new Map() as TabContent,
  outputs: new Map() as TabContent,
});

const DEFAULT_FLOAT_PRECISION = 4;

const getScalarContent = (funcCall: DG.FuncCall, prop: DG.Property) => {
  const precision = prop.options.precision;

  const scalarValue = funcCall.outputs[prop.name];
  const formattedScalarValue = prop.propertyType === DG.TYPE.FLOAT && scalarValue ?
    precision ? scalarValue.toPrecision(precision): scalarValue.toFixed(DEFAULT_FLOAT_PRECISION):
    scalarValue;
  const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

  return [scalarValue, formattedScalarValue, units] as const;
};

const tabToProperties = (fc: DG.FuncCall) => {
  const func = fc.func;
  const map = getEmptyTabToProperties();

  const processDf = (dfProp: DG.Property, isOutput: boolean) => {
    const dfViewers = Utils.getPropViewers(dfProp).config;
    if (dfViewers.length === 0) return;

    dfViewers.forEach((dfViewer) => {
      const dfBlockTitle = dfViewer.title ?? dfProp.options['caption'] ?? dfProp.name ?? ' '
      const dfNameWithViewer = `${dfBlockTitle} / ${dfViewer['type']}`;

      const tabLabel = dfProp.category === 'Misc' ?
        dfNameWithViewer: `${dfProp.category}: ${dfNameWithViewer}`;

      const name = dfProp.name;
      let df = isOutput ? fc.outputs[name] : fc.inputs[name];
      if (df)
        df = Vue.markRaw(df);
      if (isOutput)
        map.outputs.set(tabLabel, {type: 'dataframe', name, df, config: dfViewer});
      else
        map.inputs.set(tabLabel, {type: 'dataframe', name, df, config: dfViewer});
    });
    return;
  };

  func.inputs
    .forEach((inputProp) => {
      if (inputProp.propertyType === DG.TYPE.DATA_FRAME) processDf(inputProp, false);
    });

  func.outputs
    .forEach((outputProp) => {
      if (outputProp.propertyType === DG.TYPE.DATA_FRAME) {
        processDf(outputProp, true);
        return;
      }

      const category = outputProp.category === 'Misc' ? 'Output': outputProp.category;

      const categoryProps = map.outputs.get(category);
      const [rawValue, formattedValue, units] = getScalarContent(fc, outputProp);
      const scalarProp = {name: outputProp.caption || outputProp.name, rawValue, formattedValue, units}
      if (categoryProps && categoryProps.type === 'scalars')
        categoryProps.scalarsData.push(scalarProp);
      else
        map.outputs.set(category, {type: 'scalars', scalarsData: [scalarProp]});
    });

  return map;
};

const STORE_NAME = 'RFV2Layouts';

interface RFVSchema extends DBSchema {
  [STORE_NAME]: {
    key: string;
    value: PanelsState;
  };
}

export const RichFunctionView = Vue.defineComponent({
  name: 'RichFunctionView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    uuid: {
      type: String,
    },
    callState: {
      type: Object as Vue.PropType<FuncCallStateInfo>,
    },
    callMeta: {
      type: Object as Vue.PropType<Record<string, BehaviorSubject<any>>>,
    },
    validationStates: {
      type: Object as Vue.PropType<Record<string, ValidationResult>>,
    },
    consistencyStates: {
      type: Object as Vue.PropType<Record<string, ConsistencyInfo>>,
    },
    menuActions: {
      type: Object as Vue.PropType<Record<string, ViewAction[]>>,
    },
    buttonActions: {
      type: Object as Vue.PropType<ViewAction[]>,
    },
    isTreeLocked: {
      type: Boolean,
      default: false,
    },
    isReadonly: {
      type: Boolean,
      default: false,
    },
    localValidation: {
      type: Boolean,
      default: false,
    },
    historyEnabled: {
      type: Boolean,
      default: false,
    },
    showRunButton: {
      type: Boolean,
      default: true,
    },
    showStepNavigation: {
      type: Boolean,
      default: false,
    },
    skipInit: {
      type: Boolean,
      dafault: true,
    },
    viewersHook: {
      type: Function as Vue.PropType<ViewersHook>,
    },
    view: {
      type: DG.ViewBase,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (_call: DG.FuncCall) => true,
    'runClicked': () => true,
    'nextClicked': () => true,
    'actionRequested': (_actionUuid: string) => true,
    'consistencyReset': (_ioName: string) => true,
    'formValidationChanged': (_isValid: boolean) => true,
    'formInputChanged': (_a: DG.EventData<DG.InputArgs>) => true,
    'formReplaced': (_a: DG.InputForm | undefined) => true,
  },
  methods: {
    savePersonalState: () => {},
    loadPersonalLayout: () => {},
  },
  setup(props, {emit, expose}) {
    Vue.onRenderTriggered((event) => {
      console.log('RichFunctionView onRenderTriggered', event);
    });

    const {layoutDatabase} = useLayoutDb<RFVSchema>();

    const currentCall = Vue.computed(() => Vue.markRaw(props.funcCall));
    const currentView = Vue.computed(() => Vue.markRaw(props.view));
    const currentUuid = Vue.computed(() => props.uuid);
    const isFormValid = Vue.ref(false);

    const callMeta = Vue.computed(() => props.callMeta);

    const isOutputOutdated = Vue.computed(() => props.callState?.isOutputOutdated);
    const isRunning = Vue.computed(() => props.callState?.isRunning);
    const isRunnable = Vue.computed(() => props.localValidation ? isFormValid.value : props.callState?.isRunnable);
    const isReadonly = Vue.computed(() => props.isReadonly);

    const validationState = Vue.computed(() => props.validationStates);
    const consistencyState = Vue.computed(() => props.consistencyStates);

    const menuActions = Vue.computed(() => props.menuActions);
    const buttonActions = Vue.computed(() => props.buttonActions);

    const tabsData = Vue.shallowRef<RenderStateItem[]>([]);
    const tabLabels = Vue.shallowRef<string[]>([])
    const visibleTabLabels = Vue.shallowRef([] as string[]);

    const viewerTabsCount = Vue.ref<number>(0);
    const hasContextHelp = Vue.ref(false);
    const isSAenabled = Vue.ref(false);
    const isReportEnabled = Vue.ref(false);
    const isFittingEnabled = Vue.ref(false);

    const isLocked = Vue.ref(false);

    const formHidden = Vue.ref(false);
    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);

    const historyRef = Vue.shallowRef(null as InstanceType<typeof History> | null);
    const helpRef = Vue.shallowRef(null as InstanceType<typeof MarkDown> | null);
    const formRef = Vue.shallowRef(null as HTMLElement | null);
    const dockRef = Vue.shallowRef(null as InstanceType<typeof DockManager> | null);

    const layoutLoaded = Vue.ref(false);
    const dockInited = Vue.ref(false);

    const helpText = Vue.ref(null as null | string);

    const {setViewerRef} = useViewersHook(
      Vue.toRef(props, 'viewersHook'),
      Vue.toRef(props, 'callMeta'),
      currentCall,
    );

    Vue.watch(currentCall, (call) => {
      const tabToPropertiesMap = tabToProperties(call);
      tabLabels.value = [
        ...tabToPropertiesMap.inputs.keys(),
        ...tabToPropertiesMap.outputs.keys(),
      ];

      viewerTabsCount.value = [
        ...call.inputParams.values(),
        ...call.outputParams.values(),
      ].filter((param) => param.property.propertyType === DG.TYPE.DATA_FRAME)?.length;

      hasContextHelp.value = Utils.hasContextHelp(call.func);
      const features = Utils.getFeatures(call.func)
      isSAenabled.value =  Utils.getFeature(features, 'sens-analysis', false);
      isReportEnabled.value = Utils.getFeature(features, 'export', true);
      isFittingEnabled.value = Utils.getFeature(features, 'fitting', false);
    }, {immediate: true});

    Vue.watch([currentCall, isOutputOutdated, visibleTabLabels], ([call, isOutputOutdated, visibleTabLabels], [prevCall]) => {
      if (isOutputOutdated && prevCall === call)
        return;
      const tabToPropertiesMap = tabToProperties(call);

      tabsData.value = visibleTabLabels.map((tabLabel) =>
        ({
          tabLabel,
          tabContent: tabToPropertiesMap.inputs.get(tabLabel) ?? tabToPropertiesMap.outputs.get(tabLabel)!,
          isInput: !!tabToPropertiesMap.inputs.has(tabLabel)
        }));
    }, {immediate: true});

    const run = async () => {
      emit('runClicked');
    };

    const next = async () => {
      emit('nextClicked');
    };

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const tabIdx = visibleTabLabels.value.findIndex((label) => label === el.getAttribute('dock-spawn-title'));
      if (tabIdx >= 0) {
        visibleTabLabels.value.splice(tabIdx, 1);
        visibleTabLabels.value = [...visibleTabLabels.value];
      }
    };

    const personalPanelsStorage = (call: DG.FuncCall) => `${call.func.nqName}_personal_state`;
    const defaultPanelsStorage = (call: DG.FuncCall) => `${call.func.nqName}_default_state`;

    const getCurrentState = () => {
      const layout = dockRef?.value?.getLayout();
      if (!layout) return null;

      return {
        historyHidden: historyHidden.value,
        helpHidden: helpHidden.value,
        formHidden: formHidden.value,
        visibleTabLabels: visibleTabLabels.value,
        layout,
      };
    };

    const getSavedPersonalState = async (call: DG.FuncCall): Promise<PanelsState | null> => {
      const item = await layoutDatabase.value?.get(STORE_NAME, personalPanelsStorage(call));

      return item ?? null;
    };

    const getSavedDefaultState = async (call: DG.FuncCall): Promise<PanelsState | null> => {
      const item = await layoutDatabase.value?.get(STORE_NAME, defaultPanelsStorage(call));

      return item ?? null;
    };

    const saveDefaultState = async (call: DG.FuncCall = currentCall.value) => {
      if (!dockInited.value) return;

      const state = getCurrentState();
      if (state) await layoutDatabase.value?.put(STORE_NAME, state, defaultPanelsStorage(call));
    };

    const savePersonalState = async (call: DG.FuncCall = currentCall.value) => {
      if (!dockInited.value) return;

      const state = getCurrentState();
      if (state) await layoutDatabase.value?.put(STORE_NAME, state, personalPanelsStorage(call));
    };

    const removeSavedPersonalState = async () => {
      await layoutDatabase.value?.delete(STORE_NAME, personalPanelsStorage(currentCall.value));

      await loadDefaultLayout();
    };

    let intelligentLayout = true;

    const loadPersonalLayout = async () => {
      const personalState = await getSavedPersonalState(currentCall.value);

      if (!personalState)
        layoutLoaded.value = true;

      if (!dockRef.value || !personalState || !dockInited.value) return;

      intelligentLayout = false;

      historyHidden.value = personalState.historyHidden;
      helpHidden.value = personalState.helpHidden;
      formHidden.value = personalState.formHidden;
      visibleTabLabels.value = personalState.visibleTabLabels;

      await Vue.nextTick();

      await dockRef.value.useLayout(personalState.layout);

      intelligentLayout = true;
      layoutLoaded.value = true;
    };

    const loadDefaultLayout = async () => {
      const defaultState = await getSavedDefaultState(currentCall.value);
      if (!dockRef.value || !defaultState || !dockInited.value) return;

      intelligentLayout = false;

      historyHidden.value = defaultState.historyHidden;
      helpHidden.value = defaultState.helpHidden;
      formHidden.value = defaultState.formHidden;
      visibleTabLabels.value = defaultState.visibleTabLabels;

      await Vue.nextTick();

      await dockRef.value.useLayout(defaultState.layout);

      intelligentLayout = true;
    };

    Vue.watch(currentCall, async (_, oldCall) => {
      Utils.getContextHelp(currentCall.value.func).then((loadedHelp) => {
        helpText.value = loadedHelp ?? null;
      });

      if (oldCall) await savePersonalState(oldCall);
      visibleTabLabels.value = [...tabLabels.value];
    }, {immediate: true});

    const triggerSaveDefault = Vue.ref(false);
    Vue.watch(triggerSaveDefault, async () => {
      await saveDefaultState();
      await loadPersonalLayout();
    }, {flush: 'post'});

    const handleDockInit = async () => {
      layoutLoaded.value = false;
      dockInited.value = true;
      triggerSaveDefault.value = !triggerSaveDefault.value;
    };

    const onValidationChanged = (ev: boolean) => {
      if (!props.localValidation)
        return;
      isFormValid.value = ev;
      emit('formValidationChanged', ev);
    };

    const getRanges = (specificRangeName: string) => {
      const ranges: Record<string, RangeDescription> = {};
      const currentMeta = callMeta.value;
      for (const inputParam of currentCall.value.inputParams.values()) {
        const meta$ = currentMeta?.[inputParam.name];
        const range: RangeDescription = {...(meta$?.value?.[specificRangeName] ?? meta$?.value?.['range'] ?? {})};
        if (range.default == null)
          range.default = inputParam.value;
        ranges[inputParam.name] = range ?? {};
      }
      return ranges;
    };

    const getTargets = () => {
      const targets: Record<string, TargetDescription> = {};
      const currentMeta = callMeta.value;
      for (const outputParam of currentCall.value.outputParams.values()) {
        const meta$ = currentMeta?.[outputParam.name];
        const target: RangeDescription = {...(meta$?.value?.['targetFitting'] ?? {})};
        if (target.default == null)
          target.default = outputParam.value;
        targets[outputParam.name] = target ?? {};
      }
      return targets;
    };

    const runSA = () => {
      const ranges = getRanges('rangeSA');
      SensitivityAnalysisView.fromEmpty(currentCall.value.func, {ranges});
    };

    const runFitting = async () => {
      if (isLocked.value)
        return;
      isLocked.value = true;
      try {
        const currentView = grok.shell.v;
        const ranges = getRanges('rangeFitting');
        const targets = getTargets();
        const view = await FittingView.fromEmpty(currentCall.value.func, {ranges, targets, acceptMode: true});
        const call = await view.acceptedFitting$.pipe(take(1)).toPromise();
        grok.shell.v = currentView;
        if (call)
          emit('update:funcCall', Vue.markRaw(call));
      } finally {
        isLocked.value = false;
      }
    };

    expose({
      savePersonalState,
      loadPersonalLayout,
    });

    Vue.onBeforeUnmount(() => {
      savePersonalState(currentCall.value);
    });

    const menuIconStyle = {width: '15px', display: 'inline-block', textAlign: 'center'};

    return () => {
      let lastCardLabel = null as string | null;
      let scalarCardCount = 0;

      const getElementToDock = () => {
        return intelligentLayout && lastCardLabel &&
        visibleTabLabels.value.includes(lastCardLabel) &&
        scalarCardCount < 3 ? lastCardLabel: null;
      };

      const getDockRatio = () => {
        if (!intelligentLayout) return null;

        return lastCardLabel && scalarCardCount < 3 ? 0.5: 0.15;
      };

      const getDockStrategy = (categoryProps: ScalarState[]) => {
        if (!intelligentLayout) return 'fill';

        if (viewerTabsCount.value === 0) return null;

        if (categoryProps.length > 3 || scalarCardCount > 3) return 'fill';

        if (!lastCardLabel) return 'down';

        return scalarCardCount < 3 ? 'right': 'down';
      };

      return (
        Vue.withDirectives(<div class='w-full h-full flex'>
          <RibbonMenu groupName='Panels' view={currentView.value}>
            <span
              onClick={() => formHidden.value = !formHidden.value}
              class={'flex justify-between w-full'}
            >
              <div> <IconFA name='sign-in' style={menuIconStyle}/> Show inputs </div>
              { !formHidden.value && <IconFA name='check'/>}
            </span>
            <span
              onClick={() => visibleTabLabels.value = [...tabLabels.value]}
              class={'flex justify-between'}
            >
              <div> <IconFA name='sign-out'
                style={menuIconStyle}/> Show output tabs </div>
              { visibleTabLabels.value.length === tabLabels.value.length && <IconFA name='check'/>}
            </span>
            { hasContextHelp.value && <span
              onClick={() => helpHidden.value = !helpHidden.value}
              class={'flex justify-between'}
            >
              <div> <IconFA name='question' style={menuIconStyle}/> Show help </div>
              { !helpHidden.value && <IconFA name='check'/>}
            </span> }
            { props.historyEnabled && <span
              onClick={() => historyHidden.value = !historyHidden.value}
              class={'flex justify-between'}
            >
              <div> <IconFA name='history' style={menuIconStyle}/> Show history </div>
              { !historyHidden.value && <IconFA name='check'/>}
            </span> }
            <span onClick={removeSavedPersonalState}>
              <IconFA name='eraser' style={{'padding-right': '3px'}}/>
              <span> Reset layout </span>
            </span>
          </RibbonMenu>
          { menuActions.value && Object.entries(menuActions.value).map(([category, actions]) =>
            <RibbonMenu groupName={category} view={currentView.value}>
              {
                actions.map((action) => Vue.withDirectives(<span onClick={() => emit('actionRequested', action.uuid)}>
                  <div> { action.icon && <IconFA name={action.icon} style={menuIconStyle}/> } { action.friendlyName ?? action.uuid } </div>
                </span>, [[tooltip, action.description]]))
              }
            </RibbonMenu>)
          }
          <RibbonPanel view={currentView.value}>
            { isReportEnabled.value && !isOutputOutdated.value && <IconFA
              name='arrow-to-bottom'
              onClick={async () => {
                const [blob] = await Utils.richFunctionViewReport(
                  'Excel',
                  currentCall.value.func,
                  currentCall.value,
                  Utils.dfToViewerMapping(currentCall.value),
                );
                DG.Utils.download(`${currentCall.value.func.nqName} - ${Utils.getStartedOrNull(currentCall.value) ?? 'Not completed'}.xlsx`, blob);
              }}
              tooltip='Generate standard report for the current step'
            />}
            { isSAenabled.value && <IconFA
              name='analytics'
              onClick={runSA}
              tooltip='Run sensitivity analysis'
            />}
            { isFittingEnabled.value && <IconFA
              name='chart-line'
              onClick={runFitting}
              tooltip='Fit inputs'
            />}
            { hasContextHelp.value && <IconFA
              name='question'
              tooltip={ helpHidden.value ? 'Open help panel' : 'Close help panel' }
              onClick={() => helpHidden.value = !helpHidden.value}
              style={{'background-color': !helpHidden.value ? 'var(--grey-1)': null}}
            /> }
            { props.historyEnabled && <IconFA
              name='history'
              tooltip='Open history panel'
              onClick={() => historyHidden.value = !historyHidden.value}
              style={{'background-color': !historyHidden.value ? 'var(--grey-1)': null}}
            /> }
          </RibbonPanel>
          <DockManager
            key={currentUuid.value}
            onPanelClosed={handlePanelClose}
            onInitFinished={handleDockInit}
            ref={dockRef}
            class={{'pseudo_hidden': !layoutLoaded.value}}
          >
            { !historyHidden.value &&
              <History
                func={currentCall.value.func}
                showActions
                showBatchActions
                isHistory
                onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
                dock-spawn-dock-type='right'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-title='History'
                dock-spawn-panel-icon='history'
                ref={historyRef}
                class='overflow-scroll h-full'
              /> }

            { !formHidden.value &&
              <div
                class='flex flex-col p-2 overflow-scroll h-full'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-title='Inputs'
                dock-spawn-panel-icon='sign-in-alt'
                ref={formRef}
              >
                {
                  Vue.withDirectives(<InputForm
                    funcCall={currentCall.value}
                    callMeta={callMeta.value}
                    validationStates={validationState.value}
                    consistencyStates={consistencyState.value}
                    onActionRequested={(actionUuid) => emit('actionRequested', actionUuid)}
                    onConsistencyReset={(ioName) => emit('consistencyReset', ioName)}
                    onFormReplaced={(ev) =>emit('formReplaced', ev)}
                    onInputChanged={(ev) => emit('formInputChanged', ev)}
                    onValidationChanged={onValidationChanged}
                    skipInit={props.skipInit}
                    isReadonly={isReadonly.value}
                  />, [[ifOverlapping, isRunning.value, 'Recalculating...']])
                }
                <div class='flex sticky bottom-0 justify-end' style={{'z-index': 1000, 'background-color': 'rgb(255,255,255,0.75)'}}>
                  {
                    buttonActions.value?.map((action) => Vue.withDirectives(
                      <Button onClick={() => emit('actionRequested', action.uuid)}>
                        { action.icon && <IconFA name={action.icon} /> }
                        { action.friendlyName ?? action.uuid }
                      </Button>
                      , [[tooltip, action.description]]))
                  }
                  {
                    props.showRunButton &&
                      <BigButton
                        isDisabled={!isRunnable.value || isRunning.value || props.isTreeLocked || props.isReadonly}
                        onClick={run}
                      >
                        { isOutputOutdated.value ? 'Run' : 'Rerun' }
                      </BigButton>
                  }
                  {
                    props.showStepNavigation && !isOutputOutdated.value &&
                      <BigButton
                        onClick={next}
                      >
                        Next
                      </BigButton>
                  }
                </div>
              </div> }

            {
              tabsData.value
                .map(({tabLabel, tabContent, isInput}) => {
                  if (tabContent?.type === 'dataframe') {
                    const options = tabContent.config;
                    return <div
                      class='flex flex-col pl-2 h-full w-full'
                      dock-spawn-title={tabLabel}
                      dock-spawn-panel-icon={isInput ? 'sign-in-alt': 'sign-out-alt'}
                    >
                      {
                        Vue.withDirectives(<Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={tabContent.df}
                          class='w-full'
                          onViewerChanged={(v) => setViewerRef(v, tabContent.name, options['type'] as string)}
                        />, [[ifOverlapping, isRunning.value, 'Recalculating...']])
                      }
                    </div>;
                  }

                  if (tabContent?.type === 'scalars') {
                    const scalarsData = tabContent.scalarsData;

                    const panel = <ScalarsPanel
                      class='h-full overflow-scroll'
                      scalarsData={scalarsData}
                      dock-spawn-panel-icon='sign-out-alt'
                      dock-spawn-title={tabLabel}
                      dock-spawn-dock-to={getElementToDock()}
                      dock-spawn-dock-type={getDockStrategy(scalarsData)}
                      dock-spawn-dock-ratio={getDockRatio()}
                    />;

                    if (scalarsData.length < 3) {
                      scalarCardCount += scalarsData.length;
                      lastCardLabel = tabLabel;
                    }

                    return Vue.withDirectives(panel, [[ifOverlapping, isRunning.value, 'Recalculating...']]);
                  }
                })
            }
            { !helpHidden.value && helpText.value ?
              <MarkDown
                markdown={helpText.value}
                dock-spawn-title='Help'
                dock-spawn-dock-type='right'
                dock-spawn-dock-ratio={0.15}
                ref={helpRef}
              /> : null
            }
          </DockManager>
        </div>, [[ifOverlapping, isLocked.value]])
      );
    };
  },
});
