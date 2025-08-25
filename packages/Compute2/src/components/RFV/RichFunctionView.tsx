import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
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
import {DiffGrok, FittingView, TargetDescription} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {dfToViewerMapping, richFunctionViewReport, SensitivityAnalysisView} from '@datagrok-libraries/compute-utils';
import {RangeDescription} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {ScalarsPanel, ScalarState} from './ScalarsPanel';
import {BehaviorSubject} from 'rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {useViewersHook} from '../../composables/use-viewers-hook';
import {ViewAction} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {startWith, take, map} from 'rxjs/operators';
import {useHelp} from '../../composables/use-help';
import {useObservable} from '@vueuse/rxjs';
import {getIvp2WebWorker, getPipelineCreator} from '@datagrok/diff-grok';

interface ScalarsState {
  type: 'scalars',
  scalarsData: ScalarState[],
}

interface DataFrameState {
  name: string,
  df: Vue.ShallowRef<DG.DataFrame>,
  type: 'dataframe',
  config: Record<string, any>,
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
  const scalarValue = funcCall.outputs[prop.name];
  let formattedScalarValue = undefined;

  if (prop.propertyType === DG.TYPE.FLOAT && scalarValue != null) {
    if (prop.options.format)
      formattedScalarValue = DG.format(scalarValue, prop.options.format);
    else if (prop.options.precision)
      formattedScalarValue = scalarValue.toPrecision(prop.options.precision);
    else
      formattedScalarValue = scalarValue.toFixed(DEFAULT_FLOAT_PRECISION);
  }
  const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

  return [scalarValue, formattedScalarValue, units] as const;
};

const tabToProperties = (fc: DG.FuncCall) => {
  const func = fc.func;
  const tabsToProps = getEmptyTabToProperties();

  const processDf = (dfProp: DG.Property, isOutput: boolean) => {
    const dfViewers = Utils.getPropViewers(dfProp).config;
    if (dfViewers.length === 0) return;

    dfViewers.forEach((dfViewer) => {
      const dfBlockTitle = dfViewer.title ?? dfProp.options['caption'] ?? dfProp.name ?? ' ';
      const dfNameWithViewer = `${dfBlockTitle} / ${dfViewer['type']}`;

      const tabLabel = dfProp.category === 'Misc' ?
        dfNameWithViewer: `${dfProp.category}: ${dfNameWithViewer}`;

      const name = dfProp.name;
      const source = isOutput ? fc.outputParams : fc.inputParams;
      const changes$ = source[name].onChanged.pipe(
        startWith(null),
        map(() => source[name].value ? Vue.markRaw(source[name].value) : null),
      );
      const df = useObservable(changes$);
      if (isOutput)
        tabsToProps.outputs.set(tabLabel, {type: 'dataframe', name, df, config: dfViewer});
      else
        tabsToProps.inputs.set(tabLabel, {type: 'dataframe', name, df, config: dfViewer});
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

      const categoryProps = tabsToProps.outputs.get(category);
      const [rawValue, formattedValue, units] = getScalarContent(fc, outputProp);
      const scalarProp = {name: outputProp.caption || outputProp.name, rawValue, formattedValue, units};
      if (categoryProps && categoryProps.type === 'scalars')
        categoryProps.scalarsData.push(scalarProp);
      else
        tabsToProps.outputs.set(category, {type: 'scalars', scalarsData: [scalarProp]});
    });

  return tabsToProps;
};

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
  setup(props, {emit}) {
    Vue.onRenderTriggered((event) => {
      console.log('RichFunctionView onRenderTriggered', event);
    });
    const {
      helpHidden,
      helpContent,
      helpLoading,
      changeHelpFunc,
    } = useHelp();

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
    const tabLabels = Vue.shallowRef<string[]>([]);
    const visibleTabLabels = Vue.shallowRef([] as string[]);

    const isSAenabled = Vue.ref(false);
    const isReportEnabled = Vue.ref(false);
    const isFittingEnabled = Vue.ref(false);

    const isLocked = Vue.ref(false);

    const formHidden = Vue.ref(false);
    const historyHidden = Vue.ref(true);

    const historyRef = Vue.shallowRef<InstanceType<typeof History> | undefined>(undefined);
    const helpRef = Vue.shallowRef<HTMLElement | undefined>(undefined);
    const formRef = Vue.shallowRef<HTMLElement | undefined>(undefined);

    ////
    // FuncCall related
    ////

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

      const features = Utils.getFeatures(call.func);
      isSAenabled.value = Utils.getFeature(features, 'sens-analysis', false);
      isReportEnabled.value = Utils.getFeature(features, 'export', true);
      isFittingEnabled.value = Utils.getFeature(features, 'fitting', false);
    }, {immediate: true});

    Vue.watch([currentCall, isOutputOutdated, visibleTabLabels], ([call, isOutputOutdated, visibleTabLabels], [prevCall, prevOutputOutdated, prevLabels]) => {
      if (prevCall === call && visibleTabLabels === prevLabels && isOutputOutdated === prevOutputOutdated)
        return;
      const tabToPropertiesMap = tabToProperties(call);

      tabsData.value = visibleTabLabels.map((tabLabel) =>
        ({
          tabLabel,
          tabContent: tabToPropertiesMap.inputs.get(tabLabel) ?? tabToPropertiesMap.outputs.get(tabLabel)!,
          isInput: !!tabToPropertiesMap.inputs.has(tabLabel),
        }));
    }, {immediate: true});

    Vue.watch(currentCall, async (call) => {
      changeHelpFunc(call?.func);
    }, {immediate: true});

    const run = async () => {
      emit('runClicked');
    };

    const next = async () => {
      emit('nextClicked');
    };

    ////
    // DockManager related
    ////

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const tabIdx = visibleTabLabels.value.findIndex((label) => label === el.getAttribute('dock-spawn-title'));
      if (tabIdx >= 0) {
        visibleTabLabels.value.splice(tabIdx, 1);
        visibleTabLabels.value = [...visibleTabLabels.value];
      }
    };

    Vue.watch(currentCall, () => {
      visibleTabLabels.value = [...tabLabels.value];
    }, {immediate: true});

    ////
    // Intergrations related
    ////

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

    const getDiffGrok = async () => {
      let diffGrok: DiffGrok | undefined = undefined;
      let hasDiffStudio = false;
      try {
        hasDiffStudio = !!(await grok.functions.find('DiffStudio:serializeEquations')); // throws error if not found
      } catch {}
      if (hasDiffStudio) {
        const script = (currentCall.value.func as any)?.language === 'ivp' ? (currentCall.value.func as DG.Script).script : undefined;
        const ivp = await grok.functions.call('DiffStudio:serializeEquations', { problem: script });
        diffGrok = {
          ivp,
          ivpWW: getIvp2WebWorker(ivp),
          pipelineCreator: getPipelineCreator(ivp),
        };
      }
      return diffGrok;
    }

    const runSA = async () => {
      const ranges = getRanges('rangeSA');
      const diffGrok = await getDiffGrok();
      SensitivityAnalysisView.fromEmpty(currentCall.value.func, {ranges, diffGrok});
    };

    const runFitting = async () => {
      if (isLocked.value)
        return;
      isLocked.value = true;
      try {
        const currentView = grok.shell.v;
        const ranges = getRanges('rangeFitting');
        const targets = getTargets();
        const diffGrok = await getDiffGrok();
        const view = await FittingView.fromEmpty(currentCall.value.func, {ranges, targets, diffGrok, acceptMode: true});
        const call = await view.acceptedFitting$.pipe(take(1)).toPromise();
        grok.shell.v = currentView;
        if (call)
          emit('update:funcCall', Vue.markRaw(call));
      } finally {
        isLocked.value = false;
      }
    };

    ////
    // render
    ////

    const menuIconStyle = {width: '15px', display: 'inline-block', textAlign: 'center'};

    return () => (
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
              style={menuIconStyle}/> Show all outputs </div>
            { visibleTabLabels.value.length === tabLabels.value.length && <IconFA name='check'/>}
          </span>
          { <span
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
        </RibbonMenu>
        { menuActions.value && Object.entries(menuActions.value).map(([category, actions]) =>
          <RibbonMenu groupName={category} view={currentView.value}>
            {
              actions.map((action) => Vue.withDirectives(<span onClick={() => emit('actionRequested', action.uuid)}>
                <div> { action.icon && <IconFA name={action.icon} style={menuIconStyle}/> } { action.friendlyName ?? action.id } </div>
              </span>, [[tooltip, action.description]]))
            }
          </RibbonMenu>)
        }
        <RibbonPanel view={currentView.value}>
          { isReportEnabled.value && !isOutputOutdated.value && <IconFA
            name='arrow-to-bottom'
            onClick={async () => {
              const [blob] = await richFunctionViewReport(
                'Excel',
                currentCall.value.func,
                currentCall.value,
                dfToViewerMapping(currentCall.value),
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
          { <IconFA
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
        <DockManager class='block h-full'
          style={{overflow: 'hidden !important'}}
          key={currentUuid.value}
          onPanelClosed={handlePanelClose}
        >
          { !historyHidden.value && props.historyEnabled &&
              <History
                func={currentCall.value.func}
                onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
                allowCompare={true}
                forceHideInputs={false}
                showIsComplete={true}
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
                    onFormReplaced={(ev) => emit('formReplaced', ev)}
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
                        { action.friendlyName ?? action.id }
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
                    key={tabLabel}
                  >
                    {
                      Vue.withDirectives(<Viewer
                        type={options['type'] as string}
                        options={options}
                        dataFrame={tabContent.df.value}
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
                    key={tabLabel}
                  />;

                  return Vue.withDirectives(panel, [[ifOverlapping, isRunning.value, 'Recalculating...']]);
                }
              })
          }
          { !helpHidden.value ?
            <div
              dock-spawn-title='Help'
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.2}
              style={{overflow: 'scroll', height: '100%', padding: '5px'}}
              ref={helpRef}
            > { Vue.withDirectives(
                <MarkDown
                  markdown={helpContent.value ?? 'Help file is not avaliable'}
                />, [[ifOverlapping, helpLoading.value]])
              }
            </div>: null
          }
        </DockManager>
      </div>, [[ifOverlapping, isLocked.value]])
    );
  },
});
