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
  useUnwrappedCallMeta,
  useDgView,
  RibbonPanelItem, RibbonMenuItem,
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {tryParseJson} from '@datagrok-libraries/utils/src/string-utils';
import {RunBrowser} from '../RunStorage/RunBrowser';
import {KeptAlive} from '../KeptAlive/KeptAlive';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FittingView, TargetDescription} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {buildDiffGrokFromFunc} from '@datagrok-libraries/compute-utils/function-views/src/fitting/diff-studio/diff-grok-meta';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {SensitivityAnalysisView} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {RangeDescription} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {ScalarsPanel} from './ScalarsPanel';
import {BehaviorSubject} from 'rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {useViewersHook} from '../../composables/use-viewers-hook';
import {take} from 'rxjs/operators';
import {useHelp} from '../../composables/use-help';
import {useHistoryPanel} from '../../composables/use-history-panel';
import {_package} from '../../package-instance';
import {applyDefaultGridFloatFormat, canUseResults, disposeViewers, getViewers, pinView as pinViewHelper, STICKY_BAR_BACKGROUND} from '../../utils';
import {saveCallToProject, DfExportEntry} from '../../project-export';
import {getToggledInputs} from './input-toggles';
import {RenderStateItem, getEmptyTabToProperties, tabToProperties} from './rfv-tabs';
import {compareRibbonIcon} from '../RunComparison/comparison-views';


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
    // RFVApp: the form owns defaults, choices and validation; otherwise the workflow driver does
    standalone: {
      type: Boolean,
      default: false,
    },
    // run history and the save icon: 'function' lists the function's runs, 'step' only runs saved
    // for this step (emits `saveToHistory`), 'workflow' the provider's saved runs (row emits
    // `historyRunChosen`).
    history: {
      type: Object as Vue.PropType<{mode: 'function' | 'step' | 'workflow', func?: DG.Func}>,
    },
    viewersHook: {
      type: Function as Vue.PropType<ViewersHook>,
    },
    hideDefaultExport: {
      type: Boolean,
      default: false,
    },
  },
  emits: {
    'update:funcCall': (_call: DG.FuncCall) => true,
    'historyRunChosen': (_call: DG.FuncCall) => true,
    'stepRunChosen': (_call: DG.FuncCall, _readonly: boolean) => true,
    'saveToHistory': (_call: DG.FuncCall) => true,
    'compareRequested': () => true,
    'runClicked': () => true,
    'actionRequested': (_actionUuid: string, _additionalParams?: Record<string, any>) => true,
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
    const dgView = useDgView();
    const currentUuid = Vue.computed(() => props.uuid);
    const isFormValid = Vue.ref(false);

    const isOutputOutdated = Vue.computed(() => props.callState?.isOutputOutdated);
    const isRunning = Vue.computed(() => props.callState?.isRunning);
    const isRunnable = Vue.computed(() => props.standalone ? isFormValid.value : props.callState?.isRunnable);
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
    const autorun = Vue.ref(false);

    const isFittingActive = Vue.ref(false);
    const uiBlocked = Vue.computed(() => props.isBlocked || isFittingActive.value);

    const historyMode = Vue.computed(() => props.history?.mode);
    const isWorkflowHistory = Vue.computed(() => historyMode.value === 'workflow');

    const formHidden = Vue.ref(false);
    const inputToggles = Vue.shallowRef<Utils.InputToggle[]>([]);
    const shownToggles = Vue.shallowRef(new Set<string>());
    const shownTogglesKey = () => `shown_input_toggles_${currentCall.value?.func?.nqName}`;
    const switchToggle = (caption: string) => {
      const shown = new Set(shownToggles.value);
      if (!shown.delete(caption))
        shown.add(caption);
      shownToggles.value = shown;
      sessionStorage.setItem(shownTogglesKey(), JSON.stringify([...shown]));
    };
    const historyFunc = () => props.history?.func ?? currentCall.value.func;
    const {historyHidden, historyReadonly} = useHistoryPanel(() => historyFunc().nqName);
    const toggleHistory = (readonly: boolean) => {
      if (!historyHidden.value && historyReadonly.value === readonly)
        historyHidden.value = true;
      else {
        historyReadonly.value = readonly;
        historyHidden.value = false;
      }
    };

    const historyRef = Vue.shallowRef<InstanceType<typeof KeptAlive> | undefined>(undefined);

    const onHistoryRunChosen = (chosenCall: DG.FuncCall, readonly: boolean) => isWorkflowHistory.value ?
      emit('historyRunChosen', chosenCall) :
      emit('stepRunChosen', chosenCall, readonly);

    const helpRef = Vue.shallowRef<HTMLElement | undefined>(undefined);
    const formRef = Vue.shallowRef<HTMLElement | undefined>(undefined);
    const inputFormComponentRef = Vue.shallowRef<InstanceType<typeof InputForm> | undefined>(undefined);

    ////
    // FuncCall related
    ////


    const {setViewerRef} = useViewersHook(
      viewersHook,
      callMeta,
      currentCall,
    );

    // Mounted viewer per tab label; entries of closed tabs go stale (no undefined emit on
    // unmount) but are never read — collectDfExportEntries iterates visibleTabLabels only
    const liveViewers = new Map<string, DG.Viewer>();

    const hiddenByMeta = Vue.computed(() => {
      const hidden = new Set<string>();
      for (const [name, m] of Object.entries(callMetaValues)) {
        if (m && (m as any).hidden)
          hidden.add(name);
      }
      return hidden;
    });

    const toggledInputs = Vue.computed(() => getToggledInputs(
      [...currentCall.value.inputParams.values()].map((param) => param.property),
      inputToggles.value,
      shownToggles.value,
      hiddenByMeta.value,
      validationState.value,
      consistencyState.value,
    ));

    const visibleTabLabels = Vue.computed(() => tabLabels.value.filter((label) => {
      if (userClosed.value.has(label))
        return false;
      const inputContent = tabToPropertiesMap.value.inputs.get(label);
      if (inputContent && inputContent.type === 'dataframe' && hiddenByMeta.value.has(inputContent.name))
        return false;
      return true;
    }));

    // tabToProperties subscribes to param changes via useObservable, which needs an
    // active effect scope to register disposal — watcher callbacks have none, so each
    // rebuild runs in its own scope and stopping it releases the previous build's
    // subscriptions (down to the Dart-side param listeners)
    let tabsScope: Vue.EffectScope | undefined;
    const rebuildTabs = (call: DG.FuncCall) => {
      tabsScope?.stop();
      tabsScope = Vue.effectScope();
      tabToPropertiesMap.value = tabsScope.run(() => tabToProperties(call))!;
      tabLabels.value = [
        ...tabToPropertiesMap.value.inputs.keys(),
        ...tabToPropertiesMap.value.outputs.keys(),
      ];
    };
    Vue.onBeforeUnmount(() => tabsScope?.stop());

    // Per-function preferred tab, tracked separately for the input and output sides and pushed
    // to the dock as `preferredPanelTitle`. Restored only on function switch and run completion
    // (see the watcher below); a user click sets it directly so the click is respected until
    // the next restore.
    const inputKey = () => `opened_input_tab_${currentCall.value?.func?.nqName}`;
    const outputKey = () => `opened_output_tab_${currentCall.value?.func?.nqName}`;
    const preferredTab = Vue.ref<string | null>(null);

    const sideTabs = (side: 'inputs' | 'outputs') =>
      visibleTabLabels.value.filter((l) => tabToPropertiesMap.value[side].has(l));

    // Saved tab if still visible, else the default for the side (not persisted): the last output
    // tab (typically the final result) but the first input tab.
    const resolveSide = (side: 'inputs' | 'outputs', key: string) => {
      const tabs = sideTabs(side);
      const saved = sessionStorage.getItem(key);
      if (saved && tabs.includes(saved))
        return saved;
      return (side === 'outputs' ? tabs[tabs.length - 1] : tabs[0]) ?? null;
    };

    // formAsTab forces the 'Inputs' tab; otherwise input vs output side by run state.
    const resolvePreferredTab = () =>
      formAsTab.value ? 'Inputs' :
        (isOutputOutdated.value ?
          resolveSide('inputs', inputKey()) :
          resolveSide('outputs', outputKey()));

    const handleTabClicked = (title: string | null) => {
      if (!title)
        return;
      preferredTab.value = title;
      // 'Inputs' (form tab or side-panel) is not persisted: forced by formAsTab or a sticky panel.
      if (tabToPropertiesMap.value.inputs.has(title))
        sessionStorage.setItem(inputKey(), title);
      else if (tabToPropertiesMap.value.outputs.has(title))
        sessionStorage.setItem(outputKey(), title);
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
      autorun.value = props.standalone && Utils.isRunningOnInput(call.func);
      inputToggles.value = props.standalone ? [] : Utils.getInputToggles(call.func);
      shownToggles.value = new Set(tryParseJson(sessionStorage.getItem(shownTogglesKey()) ?? '[]') ?? []);
    }, {immediate: true});

    Vue.watch([currentCall, () => props.callState, visibleTabLabels], (
      [call, callState, labels], [prevCall, prevCallState, prevLabels],
    ) => {
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

      // Restore the preferred tab on function switch (incl. initial mount) and on run completion
      // (isOutputOutdated true->false). A plain visibleTabLabels change touches neither, so a
      // mid-step tab show/hide leaves focus untouched.
      const switched = prevCall !== call;
      const justRan = !switched && !!prevCallState?.isOutputOutdated && !!callState && !callState.isOutputOutdated;
      if (switched || justRan)
        preferredTab.value = resolvePreferredTab();
    }, {immediate: true});

    Vue.watch(currentCall, async (call) => {
      changeHelpFunc(call?.func);
    }, {immediate: true});

    const showRun = Vue.computed(() => !autorun.value && (isOutputOutdated.value || allowRerun.value));

    const reportHandler = async (nqName: string) => {
      try {
        await DG.Func.byName(nqName).apply({
          startDownload: true,
          funcCall: currentCall.value,
          validationState: validationState.value,
          consistencyState: consistencyState.value,
          isOutputOutdated: isOutputOutdated.value,
        });
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const exports = Vue.computed(() => {
      const activeExports: ExportItem[] = [];
      if (isReportEnabled.value && !props.hideDefaultExport) {
        const name = 'Default Excel';
        const handler = async () => {
          try {
            const viewers = await getViewers(currentCall.value, viewersHook.value, callMeta.value);
            const [blob] = await richFunctionViewReport(
              'Excel',
              currentCall.value.func,
              currentCall.value,
              viewers,
              validationState.value,
              consistencyState.value,
            ).finally(() => disposeViewers(viewers));
            DG.Utils.download(`${currentCall.value.func.nqName} - ${Utils.getStartedOrNull(currentCall.value) ?? 'Not completed'}.xlsx`, blob);
          } catch (e: any) {
            grok.shell.error(e);
          }
        };
        activeExports.push({name, handler});
      }
      activeExports.push(...customExports.value.filter((x) => x.function && x.name)
        .map((x) => ({...x, handler: () => reportHandler(x.function)})));
      activeExports.push({
        name: 'Save as project...',
        handler: () => saveCallToProject(currentCall.value, collectDfExportEntries()),
      });
      return activeExports;
    });

    // The Excel report (or a sole export) also gets a direct ribbon icon next to the menu
    const directExport = Vue.computed(() => {
      const list = exports.value;
      return list.find((e) => e.name === 'Default Excel') ?? (list.length === 1 ? list[0] : null);
    });

    // One entry per visible dataframe param, viewer configs merged across its tabs
    const collectDfExportEntries = (): DfExportEntry[] => {
      const map = tabToPropertiesMap.value;
      const entries = new Map<string, DfExportEntry>();
      for (const label of visibleTabLabels.value) {
        const isInput = map.inputs.has(label);
        const content = (isInput ? map.inputs : map.outputs).get(label);
        if (!content || content.type !== 'dataframe' || !content.df.value)
          continue;
        const key = `${isInput ? 'in' : 'out'}:${content.name}`;
        let entry = entries.get(key);
        if (!entry) {
          entry = {name: content.name, isInput, df: content.df.value, viewers: []};
          entries.set(key, entry);
        }
        const {type, ...annotated} = content.config;
        const live = liveViewers.get(label);
        const options = live ? {...annotated, ...live.getOptions().look} : annotated;
        entry.viewers.push({type: (type as string) ?? DG.VIEWER.GRID, options});
      }
      return [...entries.values()];
    };

    // On autorun the export icons stay put regardless of isOutputOutdated
    // (no ribbon rebuild/flicker); guard the click so a stale/in-flight run surfaces a shell
    // message instead of exporting partial output.
    const guardedExport = (handler: () => void) => {
      if (canUseResults(props.callState, 'exporting'))
        handler();
    };
    const exportsVisible = Vue.computed(() => autorun.value || !isOutputOutdated.value);

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

    ////
    // Intergrations related
    ////

    const onValidationChanged = (ev: boolean) => {
      if (!props.standalone)
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

    const pinView = () => pinViewHelper(dgView);

    const runSA = async () => {
      try {
        pinView();
        const ranges = getRanges('rangeSA');
        const diffGrok = await buildDiffGrokFromFunc(currentCall.value.func);
        const inputsLookup = diffGrok?.ivp?.inputsLookup ?? undefined;
        SensitivityAnalysisView.fromEmpty(currentCall.value.func,
          {ranges, diffGrok, inputsLookup, disableLookupDefault: true});
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const runFitting = async () => {
      if (isFittingActive.value)
        return;
      pinView();
      isFittingActive.value = true;
      try {
        const currentView = grok.shell.v;
        const ranges = getRanges('rangeFitting');
        const targets = getTargets();
        const diffGrok = await buildDiffGrokFromFunc(currentCall.value.func);
        const inputsLookup = diffGrok?.ivp?.inputsLookup ?? undefined;
        const view = await FittingView.fromEmpty(currentCall.value.func,
          {ranges, targets, acceptMode: true, diffGrok, inputsLookup, disableLookupDefault: true});
        const call = await view.acceptedFitting$.pipe(take(1)).toPromise();
        grok.shell.v = currentView;
        if (call)
          emit('update:funcCall', Vue.markRaw(call));
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        isFittingActive.value = false;
      }
    };

    ////
    // render
    ////

    const stepExportsItems = Vue.computed<RibbonMenuItem[]>(() => exports.value.map(({name, handler}) => ({
      text: name,
      onClick: () => guardedExport(handler),
    })));

    // the whole model or a standalone function, not a workflow step
    const isComparable = Vue.computed(() => !!props.history && historyMode.value !== 'step' && !uiBlocked.value);


    const panelsItems = Vue.computed<RibbonMenuItem[]>(() => [
      {
        text: 'Show form',
        icon: 'sign-in',
        check: !formHidden.value,
        onClick: () => formHidden.value = !formHidden.value,
      },
      {
        text: 'Show all viewers',
        icon: 'sign-out',
        check: userClosed.value.size === 0,
        onClick: () => userClosed.value = new Set(),
      },
      {
        text: 'Show help',
        icon: 'question',
        check: !helpHidden.value,
        onClick: () => helpHidden.value = !helpHidden.value,
      },
      ...(props.history ? [{
        text: 'Show history',
        icon: 'history',
        check: !historyHidden.value,
        onClick: () => toggleHistory(false),
      }] : []),
    ]);

    const ribbonItems = Vue.computed<RibbonPanelItem[]>(() => [
      ...(exportsVisible.value && !uiBlocked.value && directExport.value ? [{
        icon: 'arrow-to-bottom',
        tooltip: 'Generate report for the current step',
        onClick: () => guardedExport(directExport.value!.handler),
      }] : []),
      ...(isFittingEnabled.value && !uiBlocked.value ? [{
        icon: {path: `${_package.webRoot}files/icons/icon-chart-dots.svg`, width: 24, height: 24},
        tooltip: 'Fit inputs',
        onClick: runFitting,
      }] : []),
      ...(isSAenabled.value && !uiBlocked.value ? [{
        icon: {path: `${_package.webRoot}files/icons/icon-chart-sensitivity.svg`, width: 24, height: 24},
        tooltip: 'Run sensitivity analysis',
        onClick: runSA,
      }] : []),
      ...(props.history && (historyMode.value !== 'step' || exportsVisible.value) && !uiBlocked.value ? [{
        icon: historyMode.value === 'step' ? 'cloud-upload-alt' : 'save',
        tooltip: 'Save run',
        onClick: () => emit('saveToHistory', currentCall.value),
      }] : []),
      ...(isComparable.value ? [compareRibbonIcon(() => emit('compareRequested'))] : []),
      ...(props.history ? [{
        icon: 'history',
        tooltip: 'Open runs panel',
        active: !historyHidden.value && !historyReadonly.value,
        onClick: () => toggleHistory(false),
      }] : []),
      ...(props.history && historyMode.value === 'step' ? [{
        icon: 'lock',
        tooltip: 'Load completed run (read-only)',
        active: !historyHidden.value && historyReadonly.value,
        onClick: () => toggleHistory(true),
      }] : []),
      {
        icon: 'question',
        tooltip: helpHidden.value ? 'Open help panel' : 'Close help panel',
        active: !helpHidden.value,
        onClick: () => helpHidden.value = !helpHidden.value,
      },
    ]);

    return () => (
      Vue.withDirectives(<div class='w-full h-full flex'>
        { exportsVisible.value && !uiBlocked.value && exports.value.length > 1 &&
          <RibbonMenu groupName='Step exports' items={stepExportsItems.value}/> }
        <RibbonMenu groupName='Panels' items={panelsItems.value}/>
        <RibbonPanel items={ribbonItems.value}/>
        <DockManager class='block h-full'
          style={{overflow: 'hidden !important'}}
          onPanelClosed={handlePanelClose}
          preferredPanelTitle={preferredTab.value ?? undefined}
          onTabClicked={handleTabClicked}
          key={currentUuid.value}
        >
          { !historyHidden.value && props.history &&
            <KeptAlive
              key="__HISTORY__"
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.3}
              dock-spawn-title={historyReadonly.value ? 'Runs (read-only load)' : 'Runs'}
              dock-spawn-panel-icon='history'
              ref={historyRef}
              class='h-full'
            >
              <RunBrowser
                key={historyFunc().nqName}
                model={historyFunc()}
                stepHistory={historyMode.value === 'step'}
                openLabel={historyReadonly.value ? 'Open read-only' : 'Open'}
                onOpen={(chosenCall: DG.FuncCall) => onHistoryRunChosen(chosenCall, historyReadonly.value)}
                class='h-full'
              />
            </KeptAlive> }
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
                { [...new Set(inputToggles.value.map(({caption}) => caption))]
                  .filter((caption) => toggledInputs.value.counts.has(caption))
                  .map((caption) => {
                    const shown = shownToggles.value.has(caption);
                    const count = toggledInputs.value.counts.get(caption);
                    return <div key={caption} class='rfv2-input-toggle' onClick={() => switchToggle(caption)}>
                      <IconFA name={shown ? 'eye-slash' : 'eye'} />
                      <span>{`${shown ? 'Hide' : 'Show'} ${caption} (${count})`}</span>
                    </div>;
                  })
                }
                <InputForm
                  key={currentCall.value?.id}
                  ref={inputFormComponentRef}
                  funcCall={currentCall.value}
                  callMeta={callMeta.value}
                  validationStates={validationState.value}
                  consistencyStates={consistencyState.value}
                  onActionRequested={(actionUuid, additionalParams) =>
                    emit('actionRequested', actionUuid, additionalParams)}
                  onConsistencyReset={(ioName) => emit('consistencyReset', ioName)}
                  onFormReplaced={(ev) => emit('formReplaced', ev)}
                  onInputChanged={(ev) => emit('formInputChanged', ev)}
                  onValidationChanged={onValidationChanged}
                  skipInit={!props.standalone}
                  skipTableAutoFill={true}
                  skipLogic={!props.standalone}
                  isReadonly={isReadonly.value}
                  hiddenInputs={toggledInputs.value.hidden}
                />
                <div class='flex sticky bottom-0' style={{'z-index': 1000, 'background-color': STICKY_BAR_BACKGROUND}}>
                  { slots.navigation ?
                    slots.navigation({runLabel: runLabel.value, allowRerun: allowRerun.value}) :
                    showRun.value ?
                      <BigButton
                        isDisabled={!isRunnable.value || isRunning.value || props.isReadonly}
                        onClick={() => emit('runClicked')}
                      >
                        { isOutputOutdated.value ? runLabel.value : 'Rerun' }
                      </BigButton> :
                      (isOutputOutdated.value && !isRunnable.value && !isRunning.value) &&
                        <div
                          class='flex items-center justify-center w-full'
                          style={{color: 'var(--red-3)', fontSize: '12px', gap: '6px', padding: '6px 0'}}
                        >
                          <IconFA name='exclamation-triangle' />
                          <span>Fix invalid inputs to update results</span>
                        </div>
                  }
                </div>
              </div> }
          {
            tabsData.value
              .map(({tabLabel, tabContent, isInput}) => {
                const tabConfig = dockSpawnConfig.value[tabLabel] ?? {};
                if (tabContent?.type === 'dataframe' || tabContent?.type === 'graphics') {
                  return <div
                    class='flex flex-col pl-2 h-full w-full'
                    dock-spawn-title={tabLabel}
                    dock-spawn-panel-icon={isInput ? 'sign-in-alt': 'sign-out-alt'}
                    key={tabLabel}
                    {...tabConfig}
                  >
                    { tabContent.type === 'graphics' ?
                      (tabContent.image.value &&
                        <img
                          class={['rfv2-graphics', Utils.isSvgGraphics(tabContent.image.value) && 'rfv2-graphics-svg']}
                          src={Utils.graphicsDataUrl(tabContent.image.value)}
                        />) :
                      <Viewer
                        type={tabContent.config['type'] as string}
                        options={tabContent.config}
                        dataFrame={tabContent.df.value}
                        class='w-full'
                        onViewerChanged={(v) => {
                          if (v) liveViewers.set(tabLabel, v);
                          else liveViewers.delete(tabLabel);
                          setViewerRef(v, tabContent.name, tabContent.config['type'] as string);
                          applyDefaultGridFloatFormat(v, tabContent.config['type'] as string);
                        }}
                        onViewerDataFrameChanged={(v) =>
                          applyDefaultGridFloatFormat(v, tabContent.config['type'] as string)}
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
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.2}
              {...(dockSpawnConfig.value['Help'] ?? {})}
              dock-spawn-title='Help'
              style={{overflow: 'scroll', height: '100%', padding: '5px'}}
              key="__HELP__"
              ref={helpRef}
            > { Vue.withDirectives(
                <MarkDown
                  markdown={helpContent.value ?? 'Help file is not available'}
                />, [[ifOverlapping, helpLoading.value]])
              }
            </div>: null
          }
        </DockManager>
      </div>, [[
        ifOverlapping, isFittingActive.value, 'Fitting is active, either close it or use Apply fitted parameters',
      ]])
    );
  },
});
