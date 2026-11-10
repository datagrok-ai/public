import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {
  BigButton, Button, DockManager, IconFA, ifOverlapping, RibbonMenu, RibbonPanel, tooltip, useDgView,
  RibbonMenuItem, RibbonPanelItem,
} from '@datagrok-libraries/webcomponents-vue';
import {
  isFuncCallState,
  isStaticPipelineState,
  PipelineState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {RichFunctionView} from '../RFV/RichFunctionView';
import type {MetaCallInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {
  STEP_HISTORY_OPTION, asItemMetadata, defaultRunName,
} from '../RunStorage/run-storage';
import {openSaveRunDialog} from '../RunStorage/SaveRunDialog';
import {TreeNode} from './TreeNode';
import {StepsTree} from './StepsTree';
import {AugmentedStat} from './types';
import {PipelineView} from '../PipelineView/PipelineView';
import {useUrlSearchParams} from '@vueuse/core';
import {Inspector} from '../Inspector/Inspector';
import {
  busyReason, Control, disposeViewers, getViewers, pinView, pressControl, READONLY_REASON, runControl,
} from '../../utils';
import {
  findNextStep,
  findNextSubStep,
  findNodeWithPathByUuid, findPrevStep, findTreeNodeByPath,
  findTreeNodeParrent, getRelevantGlobalActions, hasInconsistencies, hasSubtreeFixableInconsistencies,
  hasSubtreeAnyInconsistencies,
  resolveChosenUuid, resolveSingleStep, SELECTED_STEP_BACKGROUND,
  isDeletable, isDuplicable,
} from './tree-utils';
import {
  ignoreTab, SCOPE_CLASS, shortcutsRibbonItem, useViewKeys,
} from '../../composables/use-keyboard';
import {applyCustomExport, reportTree, getExportSummary, reportSummary} from './tree-export';
import {useReactiveTreeDriver} from '../../composables/use-reactive-tree-driver';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {setHelpService} from '../../composables/use-help';
import {useFeedbackItems} from '../../composables/use-feedback';
import {createCompositorOverlayService} from '../../composables/use-compositor-overlay';
import {compositorOverlay} from '../../directives/compositor-overlay';
import {CustomExport, ExportCbInput, ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {BehaviorSubject} from 'rxjs';
import {entryFromSnapshot} from '../RunComparison/entry-extraction';
import {snapshotState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import {
  addRunToComparison, addToComparison, compareRibbonIcon, compareWithOriginal, snapshotName,
} from '../RunComparison/comparison-views';


export const TreeWizard = Vue.defineComponent({
  name: 'TreeWizard',
  props: {
    providerFunc: {
      type: String,
      required: true,
    },
    version: {
      type: String,
      required: false,
    },
    instanceConfig: {
      type: Object as Vue.PropType<PipelineInstanceConfig>,
      required: false,
    },
    initialRunId: {
      type: String,
      required: false,
    },
    showReturn: {
      type: Boolean,
      default: false,
    },
    modelName: {
      type: String,
      required: true,
    },
  },
  emits: {
    'return': (_result: any) => true,
  },
  setup(props, {emit}) {
    Vue.onRenderTriggered((event) => {
      console.log('TreeWizard onRenderTriggered', event);
    });

    const overlayService = createCompositorOverlayService();

    const {
      treeMutationsLocked,
      isGlobalLocked,
      treeState,
      currentMetaCallData,
      hasNotSavedEdits,
      states,
      logs,
      links,
      hiddenActions,
      result,
      //
      loadPipeline,
      loadAndReplaceNestedItem,
      savePipeline,
      saveDynamicItem,
      runStep,
      runSequence,
      consistencyReset,
      runAction,
      addStep,
      removeStep,
      moveStep,
      duplicateStep,
      changeFuncCall,
      returnResult,
      inspectLinks,
      inspectConfig,
    } = useReactiveTreeDriver(
      Vue.toRef(props, 'providerFunc'), Vue.toRef(props, 'version'), Vue.toRef(props, 'instanceConfig'), overlayService,
    );

    setHelpService();

    const chosenStepUuid = Vue.ref<string | undefined>();
    // compact mode: a workflow that is a single script gets no tree, no navigation, and one ribbon
    const singleStep = Vue.computed(() => treeState.value ? resolveSingleStep(treeState.value) : undefined);
    const currentView = useDgView();
    const searchParams = useUrlSearchParams<{id?: string, currentStep?: string}>('history');
    const modelName = Vue.computed(() => props.modelName);
    const providerFunc = Vue.computed(() => {
      const func = DG.Func.byName(props.providerFunc);
      return func ? Vue.markRaw(func) : undefined;
    });
    const showReturn = Vue.computed(() => props.showReturn);
    const feedbackItems = useFeedbackItems(providerFunc);

    ////
    // results
    ////

    let isResultPending = false;

    Vue.watch(result, (result) => {
      if (!isResultPending)
        return;
      isResultPending = false;
      emit('return', result);
    });

    ////
    // actions
    ////

    const runActionWithConfirmation = (uuid: string, additionalParams?: Record<string, any>) => {
      const calledAction = chosenStepState.value?.actions?.find((action) => action.uuid === uuid);
      const confirmationMessage = calledAction?.confirmationMessage;
      if (confirmationMessage) {
        ui.dialog(`Action confirmation`)
          .add(ui.markdown(confirmationMessage))
          .onOK(() => runAction(uuid, additionalParams))
          .show({center: true, modal: true});
      } else
        runAction(uuid, additionalParams);
    };

    const nodeName = (state: PipelineState) => {
      const title = states.descriptions[state.uuid]?.title;
      return typeof title === 'string' ? title : (state.friendlyName ?? state.configId);
    };

    const removeStepWithConfirm = (state: PipelineState) => {
      const name = nodeName(state);
      ui.dialog(`Remove confirmation`)
        .add(ui.divText(`Do you want to remove "${name}"? Its unsaved inputs and results will be lost.`))
        .onOK(() => removeStep(state.uuid))
        .show({center: true, modal: true});
    };

    const snapshotOf = (state: PipelineState) => entryFromSnapshot(
      snapshotState(state, states.calls, states.consistency),
      snapshotName(nodeName(state)), {
        metadata: nodeMetadata(state.uuid,
          state.uuid === treeState.value?.uuid ? currentMetaCallData.value : undefined),
        stepHistory: isFuncCallState(state) && !!state.enableHistory,
      });
    // the function whose history holds the node's saved runs
    const modelOf = (state: PipelineState) => isFuncCallState(state) ?
      state.funcCall?.func : (state.nqName ? DG.Func.byName(state.nqName) : undefined);
    const sourceTitle = (state: PipelineState) => isFuncCallState(state) ? state.sourceRunTitle : undefined;
    const compareNodeWithOriginal = (state: PipelineState) =>
      compareWithOriginal(snapshotOf(state), state.sourceRunId!, sourceTitle(state), modelOf(state));
    // an unedited saved model and read-only copies equal their saved runs; editable copies aren't tracked
    const savedRunOf = (state: PipelineState) => state.uuid === treeState.value?.uuid ?
      (hasNotSavedEdits.value ? undefined : currentMetaCallData.value?.id) :
      (state.isReadonly ? state.sourceRunId : undefined);
    const addNodeToComparison = (state: PipelineState) => {
      const runId = savedRunOf(state);
      if (runId) {
        addRunToComparison(runId, state.uuid === treeState.value?.uuid ? currentMetaCallData.value?.title :
          sourceTitle(state), modelOf(state));
      } else
        addToComparison(snapshotOf(state), modelOf(state));
    };
    const addModelToComparison = () => addNodeToComparison(treeState.value!);

    const runSubtreeWithConfirm = (startUuid: string, rerunWithConsistent?: boolean) => {
      const includeInfoInput = ui.input.bool('Also reset info-typed inputs', {value: false});
      ignoreTab(ui.dialog(`Update confirmation`)
        .add(ui.markdown(`Do you want to update input values to consistent ones and rerun substeps? You will lose inconsistent values.`))
        .add(includeInfoInput.root)
        .onOK(() => runSequence(startUuid, rerunWithConsistent, undefined, includeInfoInput.value ?? false)))
        .show({center: true, modal: true});
    };

    // a node's title, description and tags from its meta states, overridden by extra
    const nodeMetadata = (uuid?: string, extra?: Record<string, any>) => {
      const desc = uuid ? states.descriptions[uuid] : undefined;
      return {
        title: typeof desc?.title === 'string' ? desc.title : '',
        description: typeof desc?.description === 'string' ? desc.description : '',
        tags: Array.isArray(desc?.tags) ? desc.tags : [],
        ...extra,
      };
    };

    // the Save dialog starts from the loaded run's name, description and tags, then the node's links,
    // then the model name and date
    const askSaveRun = (state: PipelineState, loaded?: MetaCallInfo) => {
      const func = modelOf(state);
      if (!func)
        throw new Error(`Node ${state.uuid} has no model to save`);
      const {title, description, tags} = nodeMetadata(state.uuid);
      return openSaveRunDialog({
        func,
        name: defaultRunName(func, loaded?.title, title),
        description: loaded?.description || description,
        tags: loaded?.tags?.length ? loaded.tags : tags,
      });
    };

    const saveSubTreeState = async (uuid: string) => {
      const state = findNodeWithPathByUuid(uuid, treeState.value!)?.state;
      const result = state ? await askSaveRun(state) : null;
      if (!result)
        return;
      await saveDynamicItem(uuid, asItemMetadata(result));
    };

    const saveEntireModelState = async () => {
      if (!treeState.value) return;
      const loaded = currentMetaCallData.value;
      const result = await askSaveRun(treeState.value, loaded?.id ? loaded : undefined);
      if (!result)
        return;
      await savePipeline(asItemMetadata(result));
    };

    // a step opened by Back/Next takes the focus once the dock has placed its view, so the keyboard stays in the
    // workflow instead of falling back to the page
    let focusOpenedStep = false;
    const navigateTo = (uuid: string) => {
      focusOpenedStep = uuid !== chosenStepUuid.value;
      chosenStepUuid.value = uuid;
    };

    // chooses a step and shows its panel again if it was closed
    const showStep = (state: PipelineState) => {
      chosenStepUuid.value = state.uuid;
      if (isFuncCallState(state))
        rfvHidden.value = false;
      else
        pipelineViewHidden.value = false;
    };

    // the keyboard way to a tree row's plus: the workflow's steps to add
    const chooseStepToAdd = async (state: PipelineState) => {
      showStep(state);
      await Vue.nextTick();
      pipelineViewRef.value?.showStepsToAdd();
    };

    const onPanelOpened = (el: HTMLElement) => {
      if (focusOpenedStep && (el === rfvRef.value?.$el || el === pipelineViewRef.value?.$el)) {
        focusOpenedStep = false;
        el.focus({preventScroll: true});
      }
    };

    const goNextStep = () => {
      if (chosenStepUuid.value == null || treeState.value == null)
        return;
      const nextData = findNextStep(chosenStepUuid.value, treeState.value);
      if (nextData)
        navigateTo(nextData.state.uuid);
      else {
        const msg = `Now you can:\n
• Export results using top menu 'Export' options\n
• Save the run to your history clicking on a save icon`;
        grok.shell.info(msg);
      }
    };

    const goBack = () => {
      if (chosenStepUuid.value == null || treeState.value == null)
        return;
      const prevData = findPrevStep(chosenStepUuid.value, treeState.value);
      navigateTo(prevData ? prevData.state.uuid : treeState.value.uuid);
    };

    const onPipelineProceed = () => {
      if (chosenStepState.value && !isFuncCallState(chosenStepState.value)) {
        const nextStep = findNextSubStep(chosenStepState.value);
        if (nextStep)
          navigateTo(nextStep.state.uuid);
      }
    };

    const onRunChosen = (newCall: DG.FuncCall, readonly: boolean) => {
      if (isRootChoosen.value)
        loadPipeline(newCall.id);
      else if (chosenStepState.value && chosenStepUuid.value && treeState.value) {
        const parent = findTreeNodeParrent(chosenStepUuid.value, treeState.value);
        if (parent && !isFuncCallState(parent)) {
          const position = parent.steps.findIndex((step) => step.uuid === chosenStepUuid.value);
          loadAndReplaceNestedItem(parent.uuid, newCall.id, chosenStepState.value.configId, position, readonly);
        }
      }
    };

    const onFuncCallChange = (call: DG.FuncCall) => {
      if (chosenStepUuid.value)
        changeFuncCall(chosenStepUuid.value, call);
    };

    const onReturnClicked = () => {
      ui.dialog(`Action confirmation`)
        .add(ui.markdown(`Close this workflow and return its results`))
        .onOK(() => {
          isResultPending = true;
          returnResult();
        })
        .show({center: true, modal: true});
    };

    ////
    // DockManager related
    ////

    const treeHidden = Vue.ref(false);
    const inspectorHidden = Vue.ref(true);
    const rfvHidden = Vue.ref(true);
    const pipelineViewHidden = Vue.ref(true);
    const stepsTree = Vue.shallowRef(
      null as (InstanceType<typeof StepsTree> & {stepAt: (el: Element | null) => PipelineState | undefined}) | null);
    // the collapsed workflows outlive the tree, which unmounts while its panel is hidden
    const closedSteps = new Set<string>();
    const inspectorInstance = Vue.ref(null as InstanceType<typeof Inspector> | null);
    const rfvRef = Vue.ref(null as InstanceType<typeof RichFunctionView> | null);
    const pipelineViewRef = Vue.ref(null as (InstanceType<typeof PipelineView> & {showStepsToAdd: () => void}) | null);

    const handlePanelClose = (el: HTMLElement) => {
      if (el === stepsTree.value?.$el) treeHidden.value = true;
      if (el === inspectorInstance.value?.$el) inspectorHidden.value = true;
      if (el === rfvRef.value?.$el) rfvHidden.value = true;
      if (el === pipelineViewRef.value?.$el) pipelineViewHidden.value = true;
    };

    ////
    // routing/view integration
    ////

    const setViewName = (name: string = '') => {
      currentView.name = name;
    };

    const setViewPath = (path: string = '') => {
      currentView.path = path;
    };

    Vue.watch(searchParams, (params) => {
      const paramsRaw = [];
      if (params.currentStep)
        paramsRaw.push(`currentStep=${params.currentStep.replace(/ /g, '+')}`);
      if (params.id)
        paramsRaw.push(`id=${params.id}`);
      setViewPath(paramsRaw.length ? `?${paramsRaw.join('&')}`: '?');
    });

    const rootTitle = Vue.computed(() => nodeMetadata(treeState.value?.uuid).title);

    Vue.watch([currentMetaCallData, hasNotSavedEdits, rootTitle], ([metadata, hasNotSavedEdits, rootTitle]) => {
      if (hasNotSavedEdits)
        pinView(currentView);
      const {id, title} = metadata?.id ? metadata : {id: undefined, title: undefined};
      searchParams.id = id && !hasNotSavedEdits ? id : undefined;
      if (title)
        setViewName(`${modelName.value} - ${title}${hasNotSavedEdits ? ' (edited)' : ''}`);
      else if (!id && rootTitle)
        setViewName(`${modelName.value} - ${rootTitle}`);
      else
        setViewName(modelName.value);
    });


    const setCurrentStepByPath = (pendingStepPath: string, treeState: PipelineState) => {
      if (pendingStepPath) {
        const nodeWithPath = findTreeNodeByPath(
          pendingStepPath.split(' ').map((pathSegment) => Number.parseInt(pathSegment)),
          treeState,
        );
        if (nodeWithPath?.state)
          chosenStepUuid.value = nodeWithPath.state.uuid;
      }
    };

    let pendingStepPath: string | null = null;
    // Programmatic open (OpenWorkflowRun) passes the id via call.aux, deep links via the start URL
    let initialRunId = props.initialRunId;

    Vue.watch(treeState, (treeState) => {
      if (!treeState)
        return;

      if (pendingStepPath) {
        if (!singleStep.value)
          setCurrentStepByPath(pendingStepPath, treeState);
        pendingStepPath = null;
        return;
      }

      if (!initialRunId && globalThis.initialURLHandled)
        return;

      // Getting inital URL user entered with
      const startUrl = new URL(grok.shell.startUri);
      const loadingId = initialRunId ?? startUrl.searchParams.get('id');
      pendingStepPath = (initialRunId || singleStep.value) ? null : startUrl.searchParams.get('currentStep');
      initialRunId = undefined;
      globalThis.initialURLHandled = true;

      if (loadingId)
        loadPipeline(loadingId);
      else if (pendingStepPath) {
        setCurrentStepByPath(pendingStepPath, treeState);
        pendingStepPath = null;
      }
    }, {immediate: true});

    const chosenStep = Vue.computed(() => {
      if (!treeState.value || !chosenStepUuid.value)
        return null;
      return findNodeWithPathByUuid(chosenStepUuid.value, treeState.value) ?? null;
    });

    // Owns chosenStepUuid: defaults it on first load and repairs it when the
    // selected node disappears (e.g. the selected step was removed). Falls back
    // along the old positional path to the nearest surviving ancestor, so a
    // removed last child resolves to its parent instead of dangling.
    Vue.watch([treeState, chosenStepUuid], ([treeState]) => {
      if (!treeState)
        return;
      if (singleStep.value) {
        if (chosenStepUuid.value !== singleStep.value.step.uuid)
          chosenStepUuid.value = singleStep.value.step.uuid;
        return;
      }
      const fallbackPath = searchParams.currentStep ?
        searchParams.currentStep.split(' ').map((segment) => Number.parseInt(segment)) :
        undefined;
      const next = resolveChosenUuid(chosenStepUuid.value, treeState, fallbackPath);
      if (next !== chosenStepUuid.value)
        chosenStepUuid.value = next;
    }, {immediate: true});

    Vue.watch(chosenStep, (newStep) => {
      if (newStep && !singleStep.value)
        searchParams.currentStep = newStep.pathSegments.join(' ');
      else
        searchParams.currentStep = undefined;
    }, {immediate: true});

    ////
    // export
    ////

    const exportStates = () => ({
      callInfoStates: states.calls,
      validationStates: states.validations,
      consistencyStates: states.consistency,
      pipelineValidations: states.pipelineValidations,
      descriptions: states.descriptions,
    });

    const exports = Vue.computed(() => {
      if (!treeState.value || isFuncCallState(treeState.value))
        return [];
      const activeExports: CustomExport[] = [];
      if (!treeState.value.disableDefaultExport) {
        const defaultExport: CustomExport = {
          id: 'default',
          friendlyName: 'Default Excel',
          handler: () => reportTree({
            startDownload: true,
            treeState: treeState.value!,
            meta: currentMetaCallData.value,
            ...exportStates(),
            metaStates: states.meta,
            hasNotSavedEdits: hasNotSavedEdits.value,
          }),
        };
        activeExports.push(defaultExport);
      }
      activeExports.push(...(treeState.value.customExports ?? []));
      return activeExports;
    });

    const exportHandler = async (exportData: CustomExport) => {
      if (!treeState.value)
        return;
      const utils = {
        reportFuncCallExcel: async (
          fc: DG.FuncCall, uuid: string, viewersHook?: ViewersHook, metaState?: Record<string, BehaviorSubject<any>>,
        ) => {
          const viewers = await getViewers(fc, viewersHook, metaState);
          return richFunctionViewReport(
            'Excel',
            fc.func,
            fc,
            viewers,
            states.validations?.[uuid],
            states.consistency?.[uuid],
          ).finally(() => disposeViewers(viewers));
        },
        reportStateExcel: async (state: PipelineState, cb?: (input: ExportCbInput) => Promise<void>) => {
          return reportTree({
            startDownload: false,
            treeState: state,
            meta: currentMetaCallData.value,
            ...exportStates(),
            metaStates: states.meta,
            hasNotSavedEdits: hasNotSavedEdits.value,
            cb,
          });
        },
        getExportSummary: (state: PipelineState) => getExportSummary(state, exportStates()),
        reportSummaryExcel: (state: PipelineState) => reportSummary(getExportSummary(state, exportStates())),
        getFuncCallCustomExports: (fc: DG.FuncCall) => {
          return Utils.getCustomExports(fc.func).map((x) => x.name);
        },
        runFuncCallCustomExport: async (fc: DG.FuncCall, uuid: string, exportName: string) => {
          return applyCustomExport(fc, exportName, {
            startDownload: false,
            validationState: states.validations?.[uuid],
            consistencyState: states.consistency?.[uuid],
            isOutputOutdated: states.calls?.[uuid]?.isOutputOutdated,
            runError: states.calls?.[uuid]?.runError,
          });
        },
      };
      try {
        return await exportData.handler(treeState.value!, utils);
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const isTreeReportable = Vue.computed(() => {
      return Object.values(states.calls)
        .map((state) => state?.isOutputOutdated)
        .every((isOutdated) => isOutdated === false);
    });

    ////
    // chosen step state
    ////

    const chosenStepState = Vue.computed(() => chosenStep.value?.state);
    const isChosenParentReadonly = Vue.computed(() => !!chosenStepUuid.value && !!treeState.value &&
      !!findTreeNodeParrent(chosenStepUuid.value, treeState.value)?.isReadonly);

    // keyed on the nqName string so a new FuncCall is prepared only when the
    // chosen pipeline actually changes, not on every reactive state emission
    const chosenPipelineNqName = Vue.computed(() => {
      const state = chosenStepState.value;
      return state && !isFuncCallState(state) ? state.nqName : undefined;
    });
    const chosenPipelineFuncCall = Vue.computed(() => chosenPipelineNqName.value ?
      Vue.markRaw(DG.Func.byName(chosenPipelineNqName.value).prepare()) :
      undefined);

    // per-step history is opt-in via the `enableHistory` flag on a FuncCall step
    const currentStepHistoryEnabled = Vue.computed(() => {
      const s = chosenStepState.value;
      return !!s && isFuncCallState(s) && !!s.enableHistory;
    });

    const saveStepToHistory = async (fc: DG.FuncCall) => {
      const state = chosenStepState.value;
      const result = state ? await askSaveRun(state) : null;
      if (!result)
        return;
      fc.options[STEP_HISTORY_OPTION] = 'true';
      try {
        await historyUtils.saveRun(fc, {newId: true, metadata: result.metadata, target: result.target});
        grok.shell.info('Step saved');
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const rfvHistory = Vue.computed(() => {
      if (singleStep.value && providerFunc.value)
        return {mode: 'workflow' as const, func: providerFunc.value};
      if (currentStepHistoryEnabled.value && !isChosenParentReadonly.value)
        return {mode: 'step' as const};
      return undefined;
    });

    const isOutputOutdated = Vue.computed(() =>
      !!chosenStepUuid.value && !!states.calls[chosenStepUuid.value]?.isOutputOutdated);

    Vue.watch(chosenStepState, (state) => {
      if (!state)
        return;
      if (isFuncCallState(state)) rfvHidden.value = false;
      if (!isFuncCallState(state)) pipelineViewHidden.value = false;
    });

    const isRootChoosen = Vue.computed(() => {
      return (!!chosenStepState.value?.uuid) && chosenStepState.value?.uuid === treeState.value?.uuid;
    });

    // Ribbon controls key their visibility on tree presence only, not on the transient locks.
    // treeMutationsLocked/isGlobalLocked cycle on every link-propagation burst and every run,
    // which would otherwise hide/re-show the controls (flicker). treeState is buffered during
    // the global lock, so its last-good value keeps the controls' content stable while busy.
    const isTreeLoaded = Vue.computed(() => !!treeState.value);
    const treeBusy = Vue.computed(() => treeMutationsLocked.value || isGlobalLocked.value);

    // The ribbon lives in the view chrome, outside the compositor overlay that covers the
    // component body while busy, so the mutating controls stay visually unchanged during a lock
    // and surface a warning on click (the driver otherwise no-ops them with a console warning).
    const busyGuard = (action: string) => ({
      disabled: treeBusy.value,
      disabledReason: busyReason(action),
      disabledStyle: 'none' as const,
      disabledReasonMode: 'popup' as const,
    });

    const menuActions = Vue.computed(() => {
      if (!treeState.value || !chosenStepUuid.value)
        return {};
      const globalActions = getRelevantGlobalActions(treeState.value, chosenStepUuid.value);
      const currentStepActions = chosenStepState.value?.actions?.filter((action) => action.position === 'menu') ?? [];
      const actions = [...globalActions, ...currentStepActions]
        .filter((action) => !hiddenActions.value?.has(action.uuid));
      return actions.reduce((acc, action) => {
        const menuCategory = action.menuCategory ?? 'Actions';
        if (action.position === 'menu' || action.position === 'globalmenu') {
          if (acc[menuCategory])
            acc[menuCategory].push(action);
          else
            acc[menuCategory] = [action];
        }
        return acc;
      }, {} as Record<string, ViewAction[]>);
    });

    const buttonActions = Vue.computed(() => {
      return chosenStepState.value?.actions?.reduce((acc, action) => {
        if (action.position === 'buttons' && !hiddenActions.value?.has(action.uuid))
          acc.push(action);

        return acc;
      }, [] as ViewAction[]);
    });

    // actions on the whole workflow warn while it is busy, like its ribbon controls
    const whenIdle = (action: string, run: () => void) =>
      treeBusy.value ? grok.shell.warning(busyReason(action)) : run();
    const moveWhenIdle = (uuid: string, index: number) => whenIdle('moving a step', () => moveStep(uuid, index));
    const removeWhenIdle = (state: PipelineState) =>
      whenIdle('removing a step', () => removeStepWithConfirm(state));

    // a step's Run and Update buttons; a workflow has no Run, and its Update is the row's substeps update
    const runControlOf = (state: PipelineState): Control => {
      if (!isFuncCallState(state) || !state.funcCall)
        return {shown: false};
      const callState = states.calls[state.uuid];
      return runControl({
        isOutputOutdated: callState?.isOutputOutdated,
        allowRerun: Utils.getFeature(Utils.getFeatures(state.funcCall.func), 'rerun', false),
        autorun: false,
        isRunnable: callState?.isRunnable,
        isRunning: callState?.isRunning,
        isReadonly: state.isReadonly,
        isBlocked: treeBusy.value,
      });
    };
    const updateControlOf = (state: PipelineState): Control => ({
      shown: isFuncCallState(state) ? !!hasInconsistencies(states.consistency[state.uuid]) :
        !!hasSubtreeAnyInconsistencies(state, states.calls, states.consistency),
      disabledReason: state.isReadonly ? READONLY_REASON : treeBusy.value ? busyReason('updating') : undefined,
    });
    const updateControl = Vue.computed(() =>
      chosenStepState.value ? updateControlOf(chosenStepState.value) : {shown: false});
    const runControlOfChosen = Vue.computed(() =>
      chosenStepState.value ? runControlOf(chosenStepState.value) : {shown: false} as Control);
    const runLabel = Vue.computed(() => {
      const state = chosenStepState.value;
      const func = state && isFuncCallState(state) ? state.funcCall?.func : undefined;
      return (func ? Utils.getRunLabel(func) : null) ?? 'Run';
    });

    // Ctrl/Cmd+Enter presses Run and Ctrl/Cmd+Shift+Enter Update of the focused tree row, else of the open step;
    // Ctrl/Cmd+S saves the workflow
    useViewKeys(currentView, (e) => {
      if (!(e.ctrlKey || e.metaKey) || e.altKey)
        return;
      if (e.key === 'Enter') {
        const state = stepsTree.value?.stepAt(e.target as Element) ?? chosenStepState.value;
        if (state && (e.shiftKey ?
          pressControl(updateControlOf(state), () => runSubtreeWithConfirm(state.uuid, true)) :
          pressControl(runControlOf(state), () => runStep(state.uuid))))
          e.preventDefault();
      } else if (!e.shiftKey && e.key.toLowerCase() === 's' && isTreeLoaded.value) {
        e.preventDefault();
        whenIdle('saving', () => saveEntireModelState());
      }
    });

    ////
    // additional
    ////


    ////
    // render
    ////

    const ribbonItems = Vue.computed<RibbonPanelItem[]>(() => [
      ...(singleStep.value ? [] : [{
        icon: 'folder-tree',
        tooltip: treeHidden.value ? 'Show tree' : 'Hide tree',
        onClick: () => treeHidden.value = !treeHidden.value,
      }]),
      ...(isTreeLoaded.value && treeState.value && !singleStep.value ?
        [hasSubtreeFixableInconsistencies(treeState.value, states.calls, states.consistency) ? {
          icon: 'sync',
          tooltip: 'Update tree with consistent values',
          onClick: () => runSubtreeWithConfirm(treeState.value!.uuid, true),
          ...busyGuard('updating the tree'),
        } : {
          icon: 'forward',
          tooltip: 'Run ready steps',
          onClick: () => runSequence(treeState.value!.uuid, false),
          ...busyGuard('running the ready steps'),
        }] : []),
      ...(isTreeLoaded.value && !singleStep.value ? [{
        icon: 'save',
        tooltip: 'Save current state of model',
        onClick: () => saveEntireModelState(),
        ...busyGuard('saving'),
      }, compareRibbonIcon(addModelToComparison, busyGuard('comparing'))] : []),
      ...(isTreeLoaded.value && showReturn.value ? [{
        icon: 'check',
        tooltip: 'Confirm data',
        onClick: () => onReturnClicked(),
        ...busyGuard('confirming'),
      }] : []),
      shortcutsRibbonItem(true),
      {
        icon: 'bug',
        tooltip: inspectorHidden.value ? 'Show inspector' : 'Hide inspector',
        onClick: () => inspectorHidden.value = !inspectorHidden.value,
      },
    ]);


    const exportItems = Vue.computed<RibbonMenuItem[]>(() => exports.value.map((exportData) => ({
      text: exportData.friendlyName ?? exportData.id,
      onClick: () => exportHandler(exportData),
      ...busyGuard('exporting'),
    })));

    const actionMenuItems = (actions: ViewAction[]): RibbonMenuItem[] => actions.map((action) => ({
      text: action.friendlyName ?? action.id,
      icon: action.icon,
      tooltip: action.description,
      onClick: () => runActionWithConfirmation(action.uuid),
      ...busyGuard('running this action'),
    }));

    return () => (
      Vue.withDirectives(<div class={`w-full h-full ${SCOPE_CLASS}`}>
        <RibbonPanel items={ribbonItems.value}/>
        {isTreeLoaded.value && isTreeReportable.value && !singleStep.value &&
          <RibbonMenu groupName='Export' items={exportItems.value}/>
        }
        {feedbackItems.value.length > 0 &&
          <RibbonMenu groupName='Feedback' items={feedbackItems.value}/>
        }
        { isTreeLoaded.value && menuActions.value && Object.entries(menuActions.value).map(([category, actions]) =>
          <RibbonMenu groupName={category} items={actionMenuItems(actions)}/>)
        }
        <DockManager class='block h-full'
          style={{overflow: 'hidden !important'}}
          key={providerFunc.value?.id}
          onPanelClosed={handlePanelClose}
          onPanelOpened={onPanelOpened}
        >
          { !inspectorHidden.value &&
            <Inspector
              key="inspector"
              treeState={treeState.value}
              logs={logs.value}
              links={links.value}
              inspectLinks={inspectLinks}
              inspectConfig={inspectConfig}
              ref={inspectorInstance}
              dock-spawn-title='Inspector'
              class='h-full overflow-scroll'
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio='0.2'
            ></Inspector>
          }
          {
            treeState.value && !treeHidden.value && !singleStep.value ?
              Vue.withDirectives(<StepsTree
                key="navigation"
                dock-spawn-title='Steps'
                dock-spawn-panel-icon='folder-tree'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-z-index={2}
                treeState={treeState.value}
                closed={closedSteps}
                ref={stepsTree}
                onChoose={showStep}
                onMove={moveWhenIdle}
                onRemove={removeWhenIdle}
              >
                {({stat, isDraggable}: {stat: AugmentedStat, isDraggable: boolean}) =>
                  <TreeNode
                    stat={stat}
                    callState={states.calls[stat.data.uuid]}
                    validationStates={states.validations[stat.data.uuid]}
                    pipelineValidationState={states.pipelineValidations[stat.data.uuid]}
                    consistencyStates={states.consistency[stat.data.uuid]}
                    descriptions={states.descriptions[stat.data.uuid]}
                    style={{
                      'background-color': stat.data.uuid === chosenStepUuid.value ? SELECTED_STEP_BACKGROUND : null,
                    }}
                    isDraggable={isDraggable}
                    isDeletable={isDeletable(stat)}
                    isDuplicable={isDuplicable(stat)}
                    isReadonly={stat.data.isReadonly}
                    isSelected={stat.data.uuid === chosenStepUuid.value}
                    hasInconsistentSubsteps={
                      !!hasSubtreeAnyInconsistencies(stat.data, states.calls, states.consistency)
                    }
                    onAddNode={({itemId, position}) => addStep(stat.data.uuid, itemId, position)}
                    onChooseStepToAdd={() => chooseStepToAdd(stat.data)}
                    onRemoveNode={() => removeWhenIdle(stat.data)}
                    onDuplicateNode={() => duplicateStep(stat.data.uuid)}
                    onToggleNode={() => stat.open = !stat.open}
                    onRunSubtree={(startUuid, rerunWithConsistent) =>
                      runSubtreeWithConfirm(startUuid, rerunWithConsistent)}
                    onRunStep={(uuid) => runStep(uuid)}
                    onSaveStep={(uuid) => saveSubTreeState(uuid)}
                    onCompareWithOriginal={() => compareNodeWithOriginal(stat.data)}
                    onAddToComparison={() => addNodeToComparison(stat.data)}
                  />
                }
              </StepsTree>, [[ifOverlapping, treeMutationsLocked.value, 'Locked...']]): null
          }
          {
            !rfvHidden.value && chosenStepState.value && chosenStepUuid.value &&
            isFuncCallState(chosenStepState.value) && chosenStepState.value.funcCall &&
              <RichFunctionView
                key={chosenStepUuid.value!}
                class={{'overflow-hidden': true}}
                funcCall={chosenStepState.value.funcCall}
                uuid={chosenStepUuid.value!}
                callState={states.calls[chosenStepUuid.value]}
                callMeta={states.meta[chosenStepUuid.value]}
                viewersHook={chosenStepState.value.viewersHook}
                validationStates={states.validations[chosenStepUuid.value]}
                consistencyStates={states.consistency[chosenStepUuid.value]}
                isReadonly={chosenStepState.value.isReadonly}
                isBlocked={treeBusy.value}
                history={rfvHistory.value}
                hideDefaultExport={!!singleStep.value?.chain[0]?.disableDefaultExport}
                onUpdate:funcCall={onFuncCallChange}
                onStepRunChosen={onRunChosen}
                onHistoryRunChosen={(call) => loadPipeline(call.id)}
                onSaveToHistory={singleStep.value ? () => saveEntireModelState() : saveStepToHistory}
                onCompareRequested={addModelToComparison}
                onActionRequested={runActionWithConfirmation}
                onConsistencyReset={(ioName) => consistencyReset(chosenStepUuid.value!, ioName)}
                dock-spawn-title='Step review'
                ref={rfvRef}
              >
                {{
                  navigation: () => (
                    <>
                      { !singleStep.value &&
                        <Button onClick={goBack}>
                          Back
                        </Button>
                      }
                      {
                        buttonActions.value?.map((action) => Vue.withDirectives(
                          <Button onClick={() => runActionWithConfirmation(action.uuid)}>
                            {action.icon && <IconFA name={action.icon} />}
                            {action.friendlyName ?? action.id}
                          </Button>
                          , [[tooltip, action.description]]))
                      }
                      { runControlOfChosen.value.shown &&
                        <BigButton
                          isDisabled={!!runControlOfChosen.value.disabledReason}
                          disabledReason={runControlOfChosen.value.disabledReason}
                          onClick={() => runStep(chosenStepUuid.value!)}
                        >
                          { isOutputOutdated.value ? runLabel.value : 'Rerun' }
                        </BigButton>
                      }
                      {
                        !singleStep.value && !isOutputOutdated.value && !updateControl.value.shown &&
                          <BigButton
                            onClick={goNextStep}
                          >
                             Next
                          </BigButton>
                      }
                      { updateControl.value.shown &&
                        <BigButton
                          isDisabled={!!updateControl.value.disabledReason}
                          disabledReason={updateControl.value.disabledReason}
                          onClick={() => runSubtreeWithConfirm(chosenStepUuid.value!, true)}
                        >
                           Update
                        </BigButton>
                      }
                    </>
                  ),
                }}
              </RichFunctionView>
          }
          {
            !pipelineViewHidden.value && chosenStepUuid.value && chosenStepState.value &&
            !isFuncCallState(chosenStepState.value) &&
            <PipelineView
              funcCall={chosenPipelineFuncCall.value}
              key={chosenStepUuid.value!}
              state={chosenStepState.value}
              uuid={chosenStepUuid.value}
              isRoot={isRootChoosen.value}
              isParentReadonly={isChosenParentReadonly.value}
              buttonActions={buttonActions.value}
              body={typeof states.descriptions[chosenStepUuid.value]?.body === 'string' ?
                states.descriptions[chosenStepUuid.value]?.body as string :
                undefined}
              onActionRequested={runActionWithConfirmation}
              onCompareWithOriginal={() => compareNodeWithOriginal(chosenStepState.value!)}
              dock-spawn-title='Step sequence review'
              onUpdate:funcCall={onRunChosen}
              onAddNode={({itemId, position}) => {
                if (chosenStepUuid.value)
                  addStep(chosenStepUuid.value, itemId, position);
              }}
              ref={pipelineViewRef}
            >
              {{
                navigation: () => {
                  const isAction = isStaticPipelineState(chosenStepState.value!) && chosenStepState.value!.isActionStep;
                  return (
                    <>
                      { !isRootChoosen.value &&
                        <Button onClick={goBack}>
                          Back
                        </Button>
                      }
                      <BigButton onClick={isAction ? goNextStep : onPipelineProceed}>
                        Next
                      </BigButton>
                    </>
                  );
                },
              }}
            </PipelineView>
          }
        </DockManager>
      </div>, [[compositorOverlay, overlayService.isActive.value || !treeState.value]])
    );
  },
});
