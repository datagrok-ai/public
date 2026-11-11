import * as grok from 'datagrok-api/grok';
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
import {RunBrowser} from '../RunStorage/RunBrowser';
import {KeptAlive} from '../KeptAlive/KeptAlive';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ScalarsPanel} from './ScalarsPanel';
import {isKey, SCOPE_CLASS, shortcutsRibbonItem, TAB_CYCLE_CLASS, useViewKeys} from '../../composables/use-keyboard';
import {BehaviorSubject} from 'rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {useViewersHook} from '../../composables/use-viewers-hook';
import {useHelp} from '../../composables/use-help';
import {useHistoryPanel} from '../../composables/use-history-panel';
import {_package} from '../../package-instance';
import {applyDefaultGridFloatFormat, pressControl, runControl, STICKY_BAR_BACKGROUND} from '../../utils';
import {RenderStateItem} from './rfv-tabs';
import {compareRibbonIcon} from '../RunComparison/comparison-views';
import {useRfvFeatures} from './use-rfv-features';
import {useInputToggles} from './use-input-toggles';
import {useRfvTabs} from './use-rfv-tabs';
import {useRfvExports} from './use-rfv-exports';
import {useRfvAnalyses} from './use-rfv-analyses';


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

    const isOutputOutdated = Vue.computed(() => props.callState?.isOutputOutdated);
    const isRunning = Vue.computed(() => props.callState?.isRunning);
    const isRunnable = Vue.computed(() => props.callState?.isRunnable);
    const isReadonly = Vue.computed(() => props.isReadonly);

    const validationState = Vue.computed(() => props.validationStates);
    const consistencyState = Vue.computed(() => props.consistencyStates);

    const callMetaValues = useUnwrappedCallMeta(() => props.callMeta);

    const {
      isSAenabled, isFittingEnabled, allowRerun, runLabel, dockSpawnConfig,
      formAsTab, autorun, inputToggles,
    } = useRfvFeatures(currentCall, Vue.toRef(props, 'standalone'));

    const {isFittingActive, runSA, runFitting} = useRfvAnalyses(currentCall, callMeta, dgView,
      (call) => emit('update:funcCall', Vue.markRaw(call)));
    const uiBlocked = Vue.computed(() => props.isBlocked || isFittingActive.value);

    const historyMode = Vue.computed(() => props.history?.mode);
    const isWorkflowHistory = Vue.computed(() => historyMode.value === 'workflow');

    const formHidden = Vue.ref(false);
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

    const hiddenByMeta = Vue.computed(() => {
      const hidden = new Set<string>();
      for (const [name, m] of Object.entries(callMetaValues)) {
        if (m && (m as any).hidden)
          hidden.add(name);
      }
      return hidden;
    });

    const {shownToggles, switchToggle, toggledInputs} = useInputToggles(
      currentCall, inputToggles, hiddenByMeta, validationState, consistencyState);

    const {
      tabsData, preferredTab, userClosed, handleTabClicked, closeTab, showAllTabs, trackViewer, collectDfExportEntries,
    } = useRfvTabs(currentCall, {callState: () => props.callState, hiddenByMeta, formAsTab, isOutputOutdated});

    Vue.watch(currentCall, async (call) => {
      changeHelpFunc(call?.func);
    }, {immediate: true});

    const run = Vue.computed(() => runControl({
      isOutputOutdated: isOutputOutdated.value,
      allowRerun: allowRerun.value,
      autorun: autorun.value,
      isRunnable: isRunnable.value,
      isRunning: isRunning.value,
      isReadonly: props.isReadonly,
      isBlocked: props.isBlocked,
    }));

    const pressRun = () => pressControl(run.value, () => emit('runClicked'));

    const {exports, directExport, guardedExport, exportsVisible} = useRfvExports(currentCall, {
      hideDefaultExport: () => props.hideDefaultExport, autorun,
      callState: () => props.callState, isOutputOutdated, viewersHook, callMeta, validationState, consistencyState,
      collectDfExportEntries,
    });

    ////
    // DockManager related
    ////

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const closedLabel = el.getAttribute('dock-spawn-title');
      if (closedLabel)
        closeTab(closedLabel);
    };

    ////
    // Intergrations related
    ////

    const onValidationChanged = (ev: boolean) => {
      if (props.standalone)
        emit('formValidationChanged', ev);
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
        onClick: showAllTabs,
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
      ...(props.standalone ? [shortcutsRibbonItem(false)] : []),
    ]);

    // in a workflow the tree wizard owns the keys
    const onShortcut = (e: KeyboardEvent) => {
      if (isKey(e, ['Enter'], {mod: true})) {
        if (pressRun())
          e.preventDefault();
      } else if (isKey(e, ['s'], {mod: true}) && props.history) {
        e.preventDefault();
        if (uiBlocked.value)
          grok.shell.warning('Saving is not available while fitting is active.');
        else
          emit('saveToHistory', currentCall.value);
      }
    };
    if (props.standalone)
      useViewKeys(dgView, onShortcut);

    const inputTogglesRow = () => [...new Set(inputToggles.value.map(({caption}) => caption))]
      .filter((caption) => toggledInputs.value.counts.has(caption))
      .map((caption) => {
        const shown = shownToggles.value.has(caption);
        const count = toggledInputs.value.counts.get(caption);
        return <div key={caption} class='rfv2-input-toggle' onClick={() => switchToggle(caption)}>
          <IconFA name={shown ? 'eye-slash' : 'eye'} />
          <span>{`${shown ? 'Hide' : 'Show'} ${caption} (${count})`}</span>
        </div>;
      });

    const formPanel = () =>
      <div
        key="__FORM__"
        class={`flex flex-col p-2 overflow-scroll h-full ${TAB_CYCLE_CLASS}`}
        {...(formAsTab.value ? {} : {
          'dock-spawn-dock-type': 'left',
          'dock-spawn-dock-ratio': 0.2,
        })}
        dock-spawn-title='Inputs'
        dock-spawn-panel-icon='sign-in-alt'
        ref={formRef}
      >
        { inputTogglesRow() }
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
            slots.navigation() :
            run.value.shown ?
              <BigButton
                isDisabled={!!run.value.disabledReason}
                disabledReason={run.value.disabledReason}
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
      </div>;

    const tabPanel = ({tabLabel, tabContent, isInput}: RenderStateItem) => {
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
                trackViewer(tabLabel, v);
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
    };

    const helpPanel = () =>
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
      </div>;

    return () => (
      Vue.withDirectives(<div class={`w-full h-full flex outline-none ${SCOPE_CLASS}`} tabindex={-1}>
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
                viewersHook={historyMode.value === 'step' ? viewersHook.value : undefined}
                openLabel={historyReadonly.value ? 'Open read-only' : 'Open'}
                onOpen={(chosenCall: DG.FuncCall) => onHistoryRunChosen(chosenCall, historyReadonly.value)}
                class='h-full'
              />
            </KeptAlive> }
          { !formHidden.value && formPanel() }
          { tabsData.value.map(tabPanel) }
          { !helpHidden.value ? helpPanel() : null }
        </DockManager>
      </div>, [[
        ifOverlapping, isFittingActive.value, 'Fitting is active, either close it or use Apply fitted parameters',
      ]])
    );
  },
});
