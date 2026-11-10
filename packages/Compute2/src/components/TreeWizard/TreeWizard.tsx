import * as grok from 'datagrok-api/grok';
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
import {TreeNode} from './TreeNode';
import {StepsTree} from './StepsTree';
import {AugmentedStat} from './types';
import {PipelineView} from '../PipelineView/PipelineView';
import {Inspector} from '../Inspector/Inspector';
import {pressControl} from '../../utils';
import {
  findNextStep,
  findNextSubStep,
  findPrevStep,
  findTreeNodeParrent, hasSubtreeFixableInconsistencies,
  hasSubtreeAnyInconsistencies,
  nodeName, resolveSingleStep, SELECTED_STEP_BACKGROUND,
  isDeletable, isDuplicable,
} from './tree-utils';
import {isKey, SCOPE_CLASS, shortcutsRibbonItem, useViewKeys} from '../../composables/use-keyboard';
import {useReactiveTreeDriver} from '../../composables/use-reactive-tree-driver';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {setHelpService} from '../../composables/use-help';
import {useFeedbackItems} from '../../composables/use-feedback';
import {createCompositorOverlayService} from '../../composables/use-compositor-overlay';
import {compositorOverlay} from '../../directives/compositor-overlay';
import {compareRibbonIcon} from '../RunComparison/comparison-views';
import {confirmAction, confirmRemove, confirmReturn, confirmSubtreeUpdate} from './tree-dialogs';
import {useTreeRouting} from './use-tree-routing';
import {createNodeRuns} from './node-runs';
import {useTreeControls} from './use-tree-controls';
import {useTreeExport} from './use-tree-export';
import {useDockPanels} from './use-dock-panels';


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
    const overlayService = createCompositorOverlayService();

    const driver = useReactiveTreeDriver(
      Vue.toRef(props, 'providerFunc'), Vue.toRef(props, 'version'), Vue.toRef(props, 'instanceConfig'), overlayService,
    );
    const {
      treeMutationsLocked,
      treeState,
      states,
      logs,
      links,
      result,
      //
      loadPipeline,
      loadAndReplaceNestedItem,
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
    } = driver;

    setHelpService();

    // compact mode: a workflow that is a single script gets no tree, no navigation, and one ribbon
    const singleStep = Vue.computed(() => treeState.value ? resolveSingleStep(treeState.value) : undefined);
    const currentView = useDgView();
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

    const {chosenStepUuid, chosenStepState, chooseStep} = useTreeRouting(driver, {
      currentView, singleStep, modelName, initialRunId: props.initialRunId,
    });
    const {
      compareNodeWithOriginal, addNodeToComparison, addModelToComparison,
      saveSubTreeState, saveEntireModelState, saveStepToHistory,
    } = createNodeRuns(driver, chosenStepState);
    const {exports, exportHandler, isTreeReportable} = useTreeExport(driver);

    const stepsTree = Vue.shallowRef(
      null as (InstanceType<typeof StepsTree> & {stepAt: (el: Element | null) => PipelineState | undefined}) | null);
    // the collapsed workflows outlive the tree, which unmounts while its panel is hidden
    const closedSteps = new Set<string>();
    const inspectorInstance = Vue.ref(null as InstanceType<typeof Inspector> | null);
    const rfvRef = Vue.ref(null as InstanceType<typeof RichFunctionView> | null);
    const pipelineViewRef = Vue.ref(null as (InstanceType<typeof PipelineView> & {showStepsToAdd: () => void}) | null);
    const {
      treeHidden, inspectorHidden, rfvHidden, pipelineViewHidden,
      toggleTree, toggleInspector, handlePanelClose, navigateTo, showStep, chooseStepToAdd, onPanelOpened,
    } = useDockPanels(chosenStepUuid, chosenStepState, chooseStep, {
      stepsTree, inspector: inspectorInstance, rfv: rfvRef, pipelineView: pipelineViewRef,
    });

    const {
      isTreeLoaded, treeBusy, busyGuard, menuActions, buttonActions, whenIdle,
      runControlOf, updateControlOf, updateControl, runControlOfChosen, runLabel,
    } = useTreeControls(driver, chosenStepUuid, chosenStepState);

    ////
    // actions
    ////

    const runActionWithConfirmation = (uuid: string, additionalParams?: Record<string, any>) => confirmAction(
      chosenStepState.value?.actions?.find((action) => action.uuid === uuid)?.confirmationMessage,
      () => runAction(uuid, additionalParams));

    const removeStepWithConfirm = (state: PipelineState) =>
      confirmRemove(nodeName(states.descriptions, state), () => removeStep(state.uuid));

    const runSubtreeWithConfirm = (startUuid: string, rerunWithConsistent?: boolean) =>
      confirmSubtreeUpdate((includeInfo) => runSequence(startUuid, rerunWithConsistent, undefined, includeInfo));

    const moveWhenIdle = (uuid: string, index: number) => whenIdle('moving a step', () => moveStep(uuid, index));
    const removeWhenIdle = (state: PipelineState) =>
      whenIdle('removing a step', () => removeStepWithConfirm(state));

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

    const onReturnClicked = () => confirmReturn(() => {
      isResultPending = true;
      returnResult();
    });

    ////
    // chosen step state
    ////

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

    const rfvHistory = Vue.computed(() => {
      if (singleStep.value && providerFunc.value)
        return {mode: 'workflow' as const, func: providerFunc.value};
      if (currentStepHistoryEnabled.value && !isChosenParentReadonly.value)
        return {mode: 'step' as const};
      return undefined;
    });

    const isOutputOutdated = Vue.computed(() =>
      !!chosenStepUuid.value && !!states.calls[chosenStepUuid.value]?.isOutputOutdated);

    const isRootChoosen = Vue.computed(() => {
      return (!!chosenStepState.value?.uuid) && chosenStepState.value?.uuid === treeState.value?.uuid;
    });

    // Ctrl/Cmd+Enter presses Run and Ctrl/Cmd+Shift+Enter Update of the focused tree row, else of the open step;
    // Ctrl/Cmd+S saves the workflow; Alt+PageUp/PageDown go Back and Next, except in a grid, which pages with them
    useViewKeys(currentView, (e) => {
      if (isKey(e, ['PageUp', 'PageDown'], {alt: true}) && !singleStep.value) {
        e.preventDefault();
        if (e.key === 'PageUp')
          goBack();
        else
          goNextStep();
      } else if (isKey(e, ['Enter'], {mod: true}) || isKey(e, ['Enter'], {mod: true, shift: true})) {
        const state = stepsTree.value?.stepAt(e.target as Element) ?? chosenStepState.value;
        if (state && (e.shiftKey ?
          pressControl(updateControlOf(state), () => runSubtreeWithConfirm(state.uuid, true)) :
          pressControl(runControlOf(state), () => runStep(state.uuid))))
          e.preventDefault();
      } else if (isKey(e, ['s'], {mod: true}) && isTreeLoaded.value) {
        e.preventDefault();
        whenIdle('saving', () => saveEntireModelState());
      }
    });

    ////
    // render
    ////

    const ribbonItems = Vue.computed<RibbonPanelItem[]>(() => [
      ...(singleStep.value ? [] : [{
        icon: 'folder-tree',
        tooltip: treeHidden.value ? 'Show tree' : 'Hide tree',
        onClick: toggleTree,
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
        onClick: toggleInspector,
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

    const treeNode = ({stat, isDraggable}: {stat: AugmentedStat, isDraggable: boolean}) =>
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
      />;

    const stepNavigation = () => (
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
    );

    const pipelineNavigation = () => {
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
    };

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
                {treeNode}
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
                {{navigation: stepNavigation}}
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
              {{navigation: pipelineNavigation}}
            </PipelineView>
          }
        </DockManager>
      </div>, [[compositorOverlay, overlayService.isActive.value || !treeState.value]])
    );
  },
});
