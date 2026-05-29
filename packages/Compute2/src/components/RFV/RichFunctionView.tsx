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
  useUnwrappedCallMeta,
  DEFAULT_FLOAT_FORMAT,
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FittingView, TargetDescription} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {richFunctionViewReport, SensitivityAnalysisView} from '@datagrok-libraries/compute-utils';
import {RangeDescription} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {ScalarsPanel, ScalarsSection, ScalarState} from './ScalarsPanel';
import {BehaviorSubject} from 'rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {useViewersHook} from '../../composables/use-viewers-hook';
import {startWith, take, map} from 'rxjs/operators';
import {useHelp} from '../../composables/use-help';
import {useObservable} from '@vueuse/rxjs';
import {_package} from '../../package-instance';
import {applyDefaultGridFloatFormat, getViewers} from '../../utils';


interface ScalarsState {
  type: 'scalars',
  scalarsData: ScalarState[],
  sections?: ScalarsSection[],
}

type OutputCategoryGroupSpec = ReadonlyArray<string | Record<string, OutputCategoryGroupSpec>>;
type OutputCategoryGroups = Record<string, OutputCategoryGroupSpec>;

function parseOutputCategoryGroups(func: DG.Func): OutputCategoryGroups | undefined {
  const raw = func?.options?.['outputCategoryGroups'];
  if (!raw) return undefined;
  try {
    const parsed = typeof raw === 'string' ? JSON.parse(raw) : raw;
    return (parsed && typeof parsed === 'object' && !Array.isArray(parsed)) ? parsed : undefined;
  } catch {
    return undefined;
  }
}

function applyOutputCategoryGroups(outputs: TabContent, spec: OutputCategoryGroups) {
  const consume = (items: OutputCategoryGroupSpec, depth: number, sink: ScalarsSection[]) => {
    for (const item of items) {
      if (typeof item === 'string') {
        const entry = outputs.get(item);
        if (!entry || entry.type !== 'scalars' || entry.scalarsData.length === 0) continue;
        sink.push({label: item, indent: depth, scalarsData: entry.scalarsData});
        outputs.delete(item);
      } else if (item && typeof item === 'object' && !Array.isArray(item)) {
        for (const [label, nested] of Object.entries(item)) {
          if (!Array.isArray(nested)) continue;
          const start = sink.length;
          consume(nested, depth + 1, sink);
          if (sink.length > start)
            sink.splice(start, 0, {label, indent: depth, scalarsData: []});
        }
      }
    }
  };

  for (const [superLabel, items] of Object.entries(spec)) {
    if (!Array.isArray(items)) continue;
    const sections: ScalarsSection[] = [];
    consume(items, 0, sections);
    if (sections.length === 0) continue;
    const flat = sections.flatMap((s) => s.scalarsData);
    outputs.set(superLabel, {type: 'scalars', scalarsData: flat, sections});
  }
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

interface ExportDefinition {
  name: string,
  function: string,
}

interface ExportItem {
  name: string,
  handler: (arg?: any) => void,
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

const isEmptyScalar = (v: any) =>
  v == null || v === '' || v === DG.FLOAT_NULL || v === DG.INT_NULL;

const isEmptyDataFrame = (df: any) =>
  !df || (df instanceof DG.DataFrame && df.rowCount === 0);

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
      formattedScalarValue = DG.format(scalarValue, DEFAULT_FLOAT_FORMAT);
  } else if (typeof scalarValue === 'boolean')
    formattedScalarValue = String(scalarValue);
  const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

  return [scalarValue, formattedScalarValue, units] as const;
};

const tabToProperties = (fc: DG.FuncCall) => {
  const tabsToProps = getEmptyTabToProperties();
  const hideEmpty = !Utils.getFeature(Utils.getFeatures(fc.func), 'show-empty-outputs', false);

  const processDf = (dfProp: DG.Property, isOutput: boolean) => {
    const dfViewers = Utils.getPropViewers(dfProp).config;
    if (dfViewers.length === 0) return;
    if (hideEmpty && isOutput && isEmptyDataFrame(fc.outputs[dfProp.name])) return;

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
    if (hideEmpty && isEmptyScalar(rawValue))
      return;
    const scalarProp = {name: property.name, friendlyName: property.caption || property.name, rawValue, formattedValue, units};
    if (categoryProps && categoryProps.type === 'scalars')
      categoryProps.scalarsData.push(scalarProp);
    else
      tabsToProps.outputs.set(category, {type: 'scalars', scalarsData: [scalarProp]});
  });

  const groupsSpec = parseOutputCategoryGroups(fc.func);
  if (groupsSpec)
    applyOutputCategoryGroups(tabsToProps.outputs, groupsSpec);

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
    isBlocked: {
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
    // per-step history mode: adds a save-to-history icon (emits `saveToHistory`) and limits the
    // history panel to runs explicitly saved for this step
    stepHistory: {
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
    'saveToHistory': (_call: DG.FuncCall) => true,
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

    const viewersHook = Vue.toRef(props, 'viewersHook');
    const callMeta = Vue.toRef(props, 'callMeta');

    const currentCall = Vue.computed(() => Vue.markRaw(props.funcCall));
    const currentView = Vue.computed(() => Vue.markRaw(props.view));
    const currentUuid = Vue.computed(() => props.uuid);
    const isFormValid = Vue.ref(false);

    const isOutputOutdated = Vue.computed(() => props.callState?.isOutputOutdated);
    const isRunning = Vue.computed(() => props.callState?.isRunning);
    const isRunnable = Vue.computed(() => props.localValidation ? isFormValid.value : props.callState?.isRunnable);
    const isReadonly = Vue.computed(() => props.isReadonly);

    const validationState = Vue.computed(() => props.validationStates);
    const consistencyState = Vue.computed(() => props.consistencyStates);

    const tabsData = Vue.shallowRef<RenderStateItem[]>([]);
    const tabLabels = Vue.shallowRef<string[]>([]);
    const tabToPropertiesMap = Vue.shallowRef(getEmptyTabToProperties());
    const userClosed = Vue.shallowRef(new Set<string>());
    const callMetaValues = useUnwrappedCallMeta(() => props.callMeta);
    const dockSpawnConfig = Vue.shallowRef<Record<string, DockSpawnConfigItem>>({});
    const customExports = Vue.ref<ExportDefinition[]>([]);

    const isSAenabled = Vue.ref(false);
    const isReportEnabled = Vue.ref(false);
    const isFittingEnabled = Vue.ref(false);
    const allowRerun = Vue.ref(false);
    const runLabel = Vue.ref('Run');
    const formAsTab = Vue.ref(false);

    const isFittingActive = Vue.ref(false);
    const uiBlocked = Vue.computed(() => props.isBlocked || isFittingActive.value);

    const formHidden = Vue.ref(false);
    const inputsHidden = Vue.ref(false);
    const hasInputsHiddenOption = Vue.computed(() =>
      !!currentCall.value?.func && Utils.getInputsHidden(currentCall.value.func));
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
      viewersHook,
      callMeta,
      currentCall,
    );

    const hiddenByMeta = Vue.computed(() => {
      const hidden = new Set<string>();
      for (const [name, m] of Object.entries(callMetaValues)) {
        if (m && (m as any).hidden)
          hidden.add(name);
      }
      return hidden;
    });

    const visibleTabLabels = Vue.computed(() => tabLabels.value.filter((label) => {
      if (userClosed.value.has(label))
        return false;
      const inputContent = tabToPropertiesMap.value.inputs.get(label);
      if (inputContent && inputContent.type === 'dataframe' && hiddenByMeta.value.has(inputContent.name))
        return false;
      return true;
    }));

    const rebuildTabs = (call: DG.FuncCall) => {
      tabToPropertiesMap.value = tabToProperties(call);
      tabLabels.value = [
        ...tabToPropertiesMap.value.inputs.keys(),
        ...tabToPropertiesMap.value.outputs.keys(),
      ];
    };

    Vue.watch(currentCall, (call) => {
      rebuildTabs(call);
      userClosed.value = new Set();

      const features = Utils.getFeatures(call.func);
      isSAenabled.value = Utils.getFeature(features, 'sens-analysis', false);
      isReportEnabled.value = Utils.getFeature(features, 'export', true);
      isFittingEnabled.value = Utils.getFeature(features, 'fitting', false);
      allowRerun.value = Utils.getFeature(features, 'rerun', false);
      runLabel.value = Utils.getRunLabel(call.func) ?? 'Run';
      customExports.value = Utils.getCustomExports(call.func);
      dockSpawnConfig.value = Utils.getDockSpawnConfig(call.func);
      formAsTab.value = Utils.getFormAsTab(call.func);
      inputsHidden.value = Utils.getInputsHidden(call.func);
    }, {immediate: true});

    Vue.watch([currentCall, () => props.callState, visibleTabLabels], ([call, callState, labels], [prevCall, prevCallState, prevLabels]) => {
      if (prevCall === call && prevCallState === callState && prevLabels === labels)
        return;
      // Re-run mutates outputs in place; refresh tabToPropertiesMap and tabLabels so
      // hide-empty-outputs can drop/re-add tabs whose emptiness changed.
      if (prevCall && call === prevCall && callState !== prevCallState)
        rebuildTabs(call);
      const map = tabToPropertiesMap.value;

      tabsData.value = labels.map((tabLabel) =>
        ({
          tabLabel,
          tabContent: map.inputs.get(tabLabel) ?? map.outputs.get(tabLabel)!,
          isInput: !!map.inputs.has(tabLabel),
        }));
    }, {immediate: true});

    Vue.watch(currentCall, async (call) => {
      changeHelpFunc(call?.func);
    }, {immediate: true});

    const showRun = Vue.computed(() => props.showRunButton && (isOutputOutdated.value || allowRerun.value));

    const reportHandler = async (nqName: string) => {
      await DG.Func.byName(nqName).apply({
        startDownload: true,
        funcCall: currentCall.value,
        validationState: validationState.value,
        consistencyState: consistencyState.value,
        isOutputOutdated: isOutputOutdated.value,
      });
    }

    const exports = Vue.computed(() => {
      const activeExports: ExportItem[] = [];
      if (isReportEnabled.value) {
        const name = 'Default Excel';
        const handler = async () => {
          const viewers = await getViewers(currentCall.value, viewersHook.value, callMeta.value);
          const [blob] = await richFunctionViewReport(
            'Excel',
            currentCall.value.func,
            currentCall.value,
            viewers,
          );
          DG.Utils.download(`${currentCall.value.func.nqName} - ${Utils.getStartedOrNull(currentCall.value) ?? 'Not completed'}.xlsx`, blob);
        }
        activeExports.push({name, handler});
      }
      activeExports.push(...customExports.value.filter(x => x.function && x.name).map(x => ({...x, handler: () => reportHandler(x.function)})));
      return activeExports;
    });

    ////
    // DockManager related
    ////

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const closedLabel = el.getAttribute('dock-spawn-title');
      if (closedLabel && tabLabels.value.includes(closedLabel)) {
        const next = new Set(userClosed.value);
        next.add(closedLabel);
        userClosed.value = next;
      }
    };

    let rebuildInFlight = false;
    let rebuildTimeoutId: ReturnType<typeof setTimeout> | undefined;
    const clearRebuildFlag = () => {
      rebuildInFlight = false;
      if (rebuildTimeoutId !== undefined) {
        clearTimeout(rebuildTimeoutId);
        rebuildTimeoutId = undefined;
      }
    };
    Vue.watch(tabLabels, () => {
      rebuildInFlight = true;
      if (rebuildTimeoutId !== undefined) clearTimeout(rebuildTimeoutId);
      rebuildTimeoutId = setTimeout(clearRebuildFlag, 50);
    });
    Vue.onUnmounted(clearRebuildFlag);

    // 'Inputs' is the form side-panel unless formAsTab is on. Don't persist or
    // restore it as an active tab — it's a sticky panel, not a tab the user switches to.
    const isInputsSidePanel = (n: string | null) => n === 'Inputs' && !formAsTab.value;

    const handlePanelChanged = (name: string | null, oldName: string | null) => {
      // Restore on initial mount OR when an inflight rebuild auto-focused away from
      // the user's saved tab — push the saved tab back if it's still visible.
      if (oldName == null || rebuildInFlight) {
        const savedName = sessionStorage.getItem(`opened_tab_${currentCall.value?.func?.nqName}`);
        if (savedName && visibleTabLabels.value.includes(savedName) && !isInputsSidePanel(savedName))
          setTimeout(() => dockSpawnRef.value?.setActivePanel(savedName));
      }
      if (name && currentCall.value && !rebuildInFlight && !isInputsSidePanel(name)) {
        sessionStorage.setItem(`opened_tab_${currentCall.value.func?.nqName}`, name);
        clearRebuildFlag();
      }
    };

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
      if (isFittingActive.value)
        return;
      isFittingActive.value = true;
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
        isFittingActive.value = false;
      }
    };

    ////
    // render
    ////

    const menuIconStyle = {width: '15px', display: 'inline-block', textAlign: 'center'};

    return () => (
      Vue.withDirectives(<div class='w-full h-full flex'> { !isOutputOutdated.value && !uiBlocked.value && exports.value?.length > 1 &&
        <RibbonMenu groupName='Step exports' view={currentView.value}>
          {
            exports.value.map(({ name, handler }) =>
              <span onClick={handler}>
                <div> {name} </div>
              </span>
            )
          }
        </RibbonMenu> }
        <RibbonMenu groupName='Panels' view={currentView.value}>
          <span
            onClick={() => formHidden.value = !formHidden.value}
            class={'flex justify-between w-full'}
          >
            <div> <IconFA name='sign-in' style={menuIconStyle}/> Show form </div>
            { !formHidden.value && <IconFA name='check'/>}
          </span>
          <span
            onClick={() => userClosed.value = new Set()}
            class={'flex justify-between'}
          >
            <div> <IconFA name='sign-out'
              style={menuIconStyle}/> Show all viewers </div>
            { userClosed.value.size === 0 && <IconFA name='check'/>}
          </span>
          { <span
            onClick={() => helpHidden.value = !helpHidden.value}
            class={'flex justify-between'}
          >
            <div> <IconFA name='question' style={menuIconStyle}/> Show help </div>
            { !helpHidden.value && <IconFA name='check'/>}
          </span> }
          { (props.historyEnabled || props.stepHistory) && <span
            onClick={() => historyHidden.value = !historyHidden.value}
            class={'flex justify-between'}
          >
            <div> <IconFA name='history' style={menuIconStyle}/> Show history </div>
            { !historyHidden.value && <IconFA name='check'/>}
          </span> }
        </RibbonMenu>
        <RibbonPanel view={currentView.value}>
          { !isOutputOutdated.value && !uiBlocked.value && exports.value?.length === 1 && <IconFA
            name='arrow-to-bottom'
            onClick={exports.value[0].handler}
            tooltip='Generate report for the current step'
          /> }
          { isFittingEnabled.value && !uiBlocked.value && <IconImage
            name='fitting'
            path={`${_package.webRoot}files/icons/icon-chart-dots.svg`}
            onClick={runFitting}
            tooltip='Fit inputs'
            style={{width: '24px', height: '24px'}}
          /> }
          { isSAenabled.value && !uiBlocked.value && <IconImage
            name='sa'
            path={`${_package.webRoot}files/icons/icon-chart-sensitivity.svg`}
            onClick={runSA}
            tooltip='Run sensitivity analysis'
            style={{width: '24px', height: '24px'}}
          /> }
          { props.stepHistory && !uiBlocked.value && <IconFA
            name='cloud-upload-alt'
            tooltip='Save this step to history'
            onClick={() => emit('saveToHistory', currentCall.value)}
          /> }
          { (props.historyEnabled || props.stepHistory) && <IconFA
            name='history'
            tooltip='Open history panel'
            onClick={() => historyHidden.value = !historyHidden.value}
            style={{'background-color': !historyHidden.value ? 'var(--grey-1)': null}}
          /> }
          { <IconFA
            name='question'
            tooltip={ helpHidden.value ? 'Open help panel' : 'Close help panel' }
            onClick={() => helpHidden.value = !helpHidden.value}
            style={{'background-color': !helpHidden.value ? 'var(--grey-1)': null}}
          /> }
        </RibbonPanel>
        <DockManager class='block h-full'
          style={{overflow: 'hidden !important'}}
          onPanelClosed={handlePanelClose}
          onUpdate:activePanelTitle={handlePanelChanged}
          key={currentUuid.value}
          ref={dockSpawnRef}
        >
          { !historyHidden.value && (props.historyEnabled || props.stepHistory) &&
            <History
              key="__HISTORY__"
              func={currentCall.value.func}
              savedOnly={props.stepHistory}
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
                {...(formAsTab.value ? {} : {
                  'dock-spawn-dock-type': 'left',
                  'dock-spawn-dock-ratio': 0.2,
                })}
                dock-spawn-title='Inputs'
                dock-spawn-panel-icon='sign-in-alt'
                ref={formRef}
              >
                { hasInputsHiddenOption.value &&
                  <div
                    class='flex items-center justify-start'
                    style={{cursor: 'pointer', color: 'var(--grey-4)', fontSize: '12px', gap: '6px', padding: '4px 0'}}
                    onClick={() => inputsHidden.value = !inputsHidden.value}
                  >
                    <IconFA name={inputsHidden.value ? 'chevron-down' : 'chevron-up'} />
                    <span>{inputsHidden.value ? 'Show inputs' : 'Hide inputs'}</span>
                  </div>
                }
                { !inputsHidden.value &&
                  <InputForm
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
                  /> }
                <div class='flex sticky bottom-0' style={{'z-index': 1000, 'background-color': 'rgb(255,255,255,0.75)'}}>
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
                      <Viewer
                        type={options['type'] as string}
                        options={options}
                        dataFrame={tabContent.df.value}
                        class='w-full'
                        onViewerChanged={(v) => {
                          setViewerRef(v, tabContent.name, options['type'] as string);
                          applyDefaultGridFloatFormat(v, options['type'] as string);
                        }}
                        onViewerDataFrameChanged={(v) =>
                          applyDefaultGridFloatFormat(v, options['type'] as string)}
                      />
                    }
                  </div>;
                }

                if (tabContent?.type === 'scalars') {
                  const scalarsData = tabContent.scalarsData;
                  const sections = tabContent.sections;

                  const panel = <ScalarsPanel
                    validationStates={validationState.value}
                    class='h-full overflow-scroll'
                    scalarsData={scalarsData}
                    sections={sections}
                    dock-spawn-panel-icon='sign-out-alt'
                    dock-spawn-title={tabLabel}
                    key={tabLabel}
                    {...tabConfig}
                  />;

                  return panel;
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
      </div>, [[ifOverlapping, isFittingActive.value]])
    );
  },
});
