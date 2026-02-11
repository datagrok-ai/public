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
  IconImage,
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FittingView, TargetDescription} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {dfToViewerMapping, richFunctionViewReport, SensitivityAnalysisView} from '@datagrok-libraries/compute-utils';
import {RangeDescription} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {ScalarsPanel, ScalarState} from './ScalarsPanel';
import {BehaviorSubject} from 'rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {useViewersHook} from '../../composables/use-viewers-hook';
import {startWith, take, map} from 'rxjs/operators';
import {useHelp} from '../../composables/use-help';
import {useObservable} from '@vueuse/rxjs';
import {_package} from '../../package-instance';

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

interface DockSpawnConfigItem {
  'dock-spawn-dock-type'?: 'left' | 'right' | 'up' | 'down',
  'dock-spawn-dock-to'?: string,
  'dock-spawn-dock-ratio'?: number,
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
  const isHidden = JSON.parse(prop.options.hidden || 'false');
  if (isHidden)
    return;
  const scalarValue = funcCall.outputs[prop.name];
  let formattedScalarValue = scalarValue;

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

  [...fc.inputParams.values()].forEach(({ property }) => {
    if (property.propertyType === DG.TYPE.DATA_FRAME) processDf(property, false);
  });

  [...fc.outputParams.values()].forEach(({property}) => {
    if (property.propertyType === DG.TYPE.DATA_FRAME) {
      processDf(property, true);
      return;
    }
    const category = property.category === 'Misc' ? 'Output': property.category;

    const categoryProps = tabsToProps.outputs.get(category);
    const content = getScalarContent(fc, property);
    if (!content)
      return;
    const [rawValue, formattedValue, units] = content;
    const scalarProp = {name: property.name, friendlyName: property.caption || property.name, rawValue, formattedValue, units};
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
    'actionRequested': (_actionUuid: string) => true,
    'consistencyReset': (_ioName: string) => true,
    'formValidationChanged': (_isValid: boolean) => true,
    'formInputChanged': (_a: DG.EventData<DG.InputArgs>) => true,
    'formReplaced': (_a: DG.InputForm | undefined) => true,
  },
  setup(props, {emit, slots}) {
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

    const tabsData = Vue.shallowRef<RenderStateItem[]>([]);
    const tabLabels = Vue.shallowRef<string[]>([]);
    const visibleTabLabels = Vue.shallowRef([] as string[]);
    const activePanelTitle = Vue.shallowRef<string | undefined>(undefined);
    const dockSpawnConfig = Vue.shallowRef<Record<string, DockSpawnConfigItem>>({});

    const isSAenabled = Vue.ref(false);
    const isReportEnabled = Vue.ref(false);
    const isFittingEnabled = Vue.ref(false);
    const allowRerun = Vue.ref(false);
    const runLabel = Vue.ref('Run');

    const isLocked = Vue.ref(false);

    const formHidden = Vue.ref(false);
    const historyHidden = Vue.ref(true);

    const historyRef = Vue.shallowRef<InstanceType<typeof History> | undefined>(undefined);
    const helpRef = Vue.shallowRef<HTMLElement | undefined>(undefined);
    const formRef = Vue.shallowRef<HTMLElement | undefined>(undefined);
    const dockSpawnRef = Vue.shallowRef<InstanceType<typeof DockManager> | undefined>(undefined);
    const inputFormComponentRef = Vue.shallowRef<InstanceType<typeof InputForm> | undefined>(undefined);

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
      allowRerun.value = Utils.getFeature(features, 'rerun', false);
      runLabel.value = Utils.getRunLabel(call.func) ?? 'Run';
      dockSpawnConfig.value = Utils.getDockSpawnConfig(call.func);
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

    const showRun = Vue.computed(() => props.showRunButton && (isOutputOutdated.value || allowRerun.value));

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

    const handlePanelChanged = (name: string | null, oldName: string | null) => {
      if (oldName == null) {
        const savedName = sessionStorage.getItem(`opened_tab_${currentCall.value.func?.nqName}`);
        if (savedName)
          setTimeout(() => activePanelTitle.value = savedName);
      }

      if (currentCall.value)
        sessionStorage.setItem(`opened_tab_${currentCall.value.func?.nqName}`, name ?? '');
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
          /> }
          { isFittingEnabled.value && <IconImage
            name='fitting'
            path={`${_package.webRoot}files/icons/icon-chart-dots.svg`}
            onClick={runFitting}
            tooltip='Fit inputs'
            style={{width: '24px', height: '24px'}}
          /> }
          { isSAenabled.value && <IconImage
            name='sa'
            path={`${_package.webRoot}files/icons/icon-chart-sensitivity.svg`}
            onClick={runSA}
            tooltip='Run sensitivity analysis'
            style={{width: '24px', height: '24px'}}
          /> }
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
          onPanelClosed={handlePanelClose}
          onUpdate:activePanelTitle={handlePanelChanged}
          key={currentUuid.value}
          activePanelTitle={activePanelTitle.value}
          ref={dockSpawnRef}
        >
          { !historyHidden.value && props.historyEnabled &&
            <History
              key="__HISTORY__"
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
                key="__FORM__"
                class='flex flex-col p-2 overflow-scroll h-full'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-title='Inputs'
                dock-spawn-panel-icon='sign-in-alt'
                ref={formRef}
              >
                {
                  Vue.withDirectives(<InputForm
                    key={currentCall.value?.id}
                    ref={inputFormComponentRef}
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
                  { slots.navigation ?
                    slots.navigation({runLabel: runLabel.value, allowRerun: allowRerun.value}) :
                    showRun.value &&
                      <BigButton
                        isDisabled={!isRunnable.value || isRunning.value || props.isReadonly}
                        onClick={() => emit('runClicked')}
                      >
                        { isOutputOutdated.value ? runLabel.value : 'Rerun' }
                      </BigButton>
                  }
                </div>
              </div> }
          {
            tabsData.value
              .map(({tabLabel, tabContent, isInput}) => {
                const tabConfig = dockSpawnConfig.value[tabLabel] ?? {};
                if (tabContent?.type === 'dataframe') {
                  const options = tabContent.config;
                  return <div
                    class='flex flex-col pl-2 h-full w-full'
                    dock-spawn-title={tabLabel}
                    dock-spawn-panel-icon={isInput ? 'sign-in-alt': 'sign-out-alt'}
                    key={tabLabel}
                    {...tabConfig}
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
                    validationStates={validationState.value}
                    class='h-full overflow-scroll'
                    scalarsData={scalarsData}
                    dock-spawn-panel-icon='sign-out-alt'
                    dock-spawn-title={tabLabel}
                    key={tabLabel}
                    {...tabConfig}
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
              key="__HELP__"
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
