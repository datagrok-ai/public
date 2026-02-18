import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BigButton, Button, DockManager, IconFA, ifOverlapping, RibbonMenu, RibbonPanel, tooltip} from '@datagrok-libraries/webcomponents-vue';
import {
  isFuncCallState, isParallelPipelineState,
  isSequentialPipelineState,
  PipelineState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {RichFunctionView} from '../RFV/RichFunctionView';
import {TreeNode} from './TreeNode';
import {Draggable, dragContext} from '@he-tree/vue';
import {AugmentedStat} from './types';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {PipelineView} from '../PipelineView/PipelineView';
import {useUrlSearchParams} from '@vueuse/core';
import {Inspector} from '../Inspector/Inspector';
import {
  findNextStep,
  findNextSubStep,
  findNodeWithPathByUuid, findPrevStep, findTreeNodeByPath,
  findTreeNodeParrent, getRelevantGlobalActions, hasInconsistencies, hasSubtreeFixableInconsistencies,
  reportTree,
} from '../../utils';
import {useReactiveTreeDriver} from '../../composables/use-reactive-tree-driver';
import {take} from 'rxjs/operators';
import {EditRunMetadataDialog} from '@datagrok-libraries/compute-utils/shared-components/src/history-dialogs';
import {PipelineInstanceConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {setHelpService} from '../../composables/use-help';
import {CustomExport, ExportCbInput} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {dfToViewerMapping, richFunctionViewReport} from '@datagrok-libraries/compute-utils';

const DEVELOPERS_GROUP = 'Developers';

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
    showReturn: {
      type: Boolean,
      default: false,
    },
    modelName: {
      type: String,
      required: true,
    },
    view: {
      type: DG.ViewBase,
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

    const {
      treeMutationsLocked,
      isGlobalLocked,
      treeState,
      currentMetaCallData,
      hasNotSavedEdits,
      states,
      logs,
      config,
      links,
      result,
      //
      loadPipeline,
      loadAndReplaceNestedPipeline,
      savePipeline,
      saveDynamicItem,
      runStep,
      runSequence,
      consistencyReset,
      runAction,
      addStep,
      removeStep,
      moveStep,
      changeFuncCall,
      returnResult,
    } = useReactiveTreeDriver(Vue.toRef(props, 'providerFunc'), Vue.toRef(props, 'version'), Vue.toRef(props, 'instanceConfig'));

    setHelpService();

    const chosenStepUuid = Vue.ref<string | undefined>();
    const currentView = Vue.computed(() => Vue.markRaw(props.view));
    const searchParams = useUrlSearchParams<{id?: string, currentStep?: string}>('history');
    const modelName = Vue.computed(() => props.modelName);
    const isTreeReady = Vue.computed(() => treeState.value && !treeMutationsLocked.value && !isGlobalLocked.value);
    const providerFunc = Vue.computed(() => {
      const func = DG.Func.byName(props.providerFunc);
      return func ? Vue.markRaw(func) : undefined;
    });
    const showReturn = Vue.computed(() => props.showReturn);
    const reportBugUrl = Vue.computed<string | undefined>(() => {
      return (providerFunc.value?.package?.settings)?.REPORT_BUG_URL;
    });
    const reqFeatureUrl = Vue.computed<string | undefined>(() => {
      return (providerFunc.value?.package?.settings)?.REQUEST_FEATURE_URL;
    });

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

    const runActionWithConfirmation = (uuid: string) => {
      const calledAction = chosenStepState.value?.actions?.find((action) => action.uuid === uuid);
      const confirmationMessage = calledAction?.confirmationMessage;
      if (confirmationMessage) {
        ui.dialog(`Action confirmation`)
          .add(ui.markdown(confirmationMessage))
          .onOK(() => runAction(uuid))
          .show({center: true, modal: true});
      } else
        runAction(uuid);
    };

    const runSubtreeWithConfirm = (startUuid: string, rerunWithConsistent?: boolean) => {
      ui.dialog(`Update confirmation`)
        .add(ui.markdown(`Do you want to update input values to consistent ones and rerun substeps? You will lose inconsistent values.`))
        .onOK(() => runSequence(startUuid, rerunWithConsistent))
        .show({center: true, modal: true});
    };

    const saveSubTreeState = (uuid: string) => {
      const chosenStepDesc = states.descriptions[uuid];
      const dialog = new EditRunMetadataDialog({
        title: typeof(chosenStepDesc?.title) === 'string' ? chosenStepDesc?.title : '',
        description: typeof(chosenStepDesc?.description) === 'string' ? chosenStepDesc?.description : '',
        tags: Array.isArray(chosenStepDesc?.tags) ? chosenStepDesc?.tags : [],
      });
      dialog.onMetadataEdit.pipe(take(1)).subscribe((editOptions) => {
        saveDynamicItem(chosenStepUuid.value!, editOptions);
      });

      dialog.show({center: true, width: 500});
    };

    const saveEntireModelState = () => {
      if (!treeState.value) return;

      const rootDesc = states.descriptions[treeState.value.uuid];
      const dialog = new EditRunMetadataDialog({...rootDesc, ...currentMetaCallData.value});
      dialog.onMetadataEdit.pipe(take(1)).subscribe((editOptions) => {
        savePipeline(editOptions);
      });

      dialog.show({center: true, width: 500});
    };

    const goNextStep = () => {
      if (chosenStepUuid.value == null || treeState.value == null)
        return;
      const nextData = findNextStep(chosenStepUuid.value, treeState.value);
      if (nextData)
        chosenStepUuid.value = nextData.state.uuid;
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
      if (prevData)
        chosenStepUuid.value = prevData.state.uuid;
      else
        chosenStepUuid.value = treeState.value.uuid;
    };

    const onPipelineProceed = () => {
      if (chosenStepState.value && !isFuncCallState(chosenStepState.value)) {
        const nextStep = findNextSubStep(chosenStepState.value);
        if (nextStep)
          chosenStepUuid.value = nextStep.state.uuid;
      }
    };

    const onPipelineFuncCallUpdate = (newCall: DG.FuncCall) => {
      if (isRootChoosen.value)
        loadPipeline(newCall.id);
      else if (chosenStepState.value && chosenStepUuid.value && treeState.value && !isFuncCallState(chosenStepState.value)) {
        const parent = findTreeNodeParrent(chosenStepUuid.value, treeState.value);
        if (parent && !isFuncCallState(parent)) {
          const position = parent.steps.findIndex((step) => step.uuid === chosenStepUuid.value);
          loadAndReplaceNestedPipeline(parent.uuid, newCall.id, chosenStepState.value.configId, position);
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
    const treeInstance = Vue.shallowRef(null as InstanceType<typeof Draggable> | null);
    const inspectorInstance = Vue.ref(null as InstanceType<typeof Inspector> | null);
    const rfvRef = Vue.ref(null as InstanceType<typeof RichFunctionView> | null);
    const pipelineViewRef = Vue.ref(null as InstanceType<typeof PipelineView> | null);

    const handlePanelClose = (el: HTMLElement) => {
      if (el === treeInstance.value?.$el) treeHidden.value = true;
      if (el === inspectorInstance.value?.$el) inspectorHidden.value = true;
      if (el === rfvRef.value?.$el) rfvHidden.value = true;
      if (el === pipelineViewRef.value?.$el) pipelineViewHidden.value = true;
    };

    ////
    // routing/view integration
    ////

    const setViewName = (name: string = '') => {
      if (props.view)
        props.view.name = name;
    };

    const setViewPath = (path: string = '') => {
      if (props.view)
        props.view.path = path;
    };

    Vue.watch(searchParams, (params) => {
      const paramsRaw = [];
      if (params.currentStep)
        paramsRaw.push(`currentStep=${params.currentStep.replace(' ', '+')}`);
      if (params.id)
        paramsRaw.push(`id=${params.id}`);
      setViewPath(paramsRaw.length ? `?${paramsRaw.join('&')}`: '?');
    });

    Vue.watch([currentMetaCallData, hasNotSavedEdits], ([metadata, hasNotSavedEdits]) => {
      if (!metadata || hasNotSavedEdits) {
        searchParams.id = undefined;
        setViewName(modelName.value);
        return;
      }

      const {id, title, started} = metadata;
      if (id) searchParams.id = id;
      if (title) setViewName(`${modelName.value} - ${title}`);
      else if (started) setViewName(`${modelName.value} - ${started}`);
      else setViewName(modelName.value);
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

    Vue.watch(treeState, (treeState) => {
      if (!treeState)
        return;

      if (pendingStepPath) {
        setCurrentStepByPath(pendingStepPath, treeState);
        pendingStepPath = null;
        return;
      }

      if (globalThis.initialURLHandled)
        return;

      // Getting inital URL user entered with
      const startUrl = new URL(grok.shell.startUri);
      globalThis.initialURLHandled = true;

      const loadingId = startUrl.searchParams.get('id');
      pendingStepPath = startUrl.searchParams.get('currentStep');
      if (loadingId)
        loadPipeline(loadingId);
      else if (pendingStepPath) {
        setCurrentStepByPath(pendingStepPath, treeState);
        pendingStepPath = null;
      }
    }, {immediate: true});

    const chosenStep = Vue.computed(() => {
      if (!treeState.value)
        return null;

      if (!chosenStepUuid.value)
        chosenStepUuid.value = treeState.value.uuid;

      const step = findNodeWithPathByUuid(chosenStepUuid.value, treeState.value);

      if (step)
        return step;

      // step with current uuid is removed from the tree, try setting the current step by the route
      const currentURL = new URL(window.location.href);
      const currentStep = currentURL.searchParams.get('currentStep');
      if (currentStep)
        setCurrentStepByPath(currentStep, treeState.value);

      return findNodeWithPathByUuid(chosenStepUuid.value, treeState.value);
    });

    Vue.watch(chosenStep, (newStep) => {
      if (newStep)
        searchParams.currentStep = newStep.pathSegments.join(' ');
      else
        searchParams.currentStep = undefined;
    }, {immediate: true});

    ////
    // export
    ////

    const exports = Vue.computed(() => {
      if (!treeState.value || isFuncCallState(treeState.value))
        return [];
      const defaultExport: CustomExport = {
        id: 'default',
        friendlyName: 'Default Excel',
        handler: () => reportTree({
          startDownload: true,
          treeState: treeState.value!,
          meta: currentMetaCallData.value,
          callInfoStates: states.calls,
          validationStates: states.validations,
          consistencyStates: states.consistency,
          descriptions: states.descriptions,
          hasNotSavedEdits: hasNotSavedEdits.value
        })};
      return [defaultExport,...(treeState.value.customExports ?? [])];
    });

    const exportHandler = async (exportData: CustomExport) => {
      if (!treeState.value)
        return
      const utils = {
        reportFuncCallExcel: async (fc: DG.FuncCall, uuid: string) => {
          return richFunctionViewReport(
            'Excel',
            fc.func,
            fc,
            dfToViewerMapping(fc),
            states.validations?.[uuid],
            states.consistency?.[uuid],
          );
        },
        reportStateExcel: async (state: PipelineState, cb?: (input: ExportCbInput) => Promise<void>) => {
          return reportTree({
            startDownload: false,
            treeState: state,
            meta: currentMetaCallData.value,
            callInfoStates: states.calls,
            validationStates: states.validations,
            consistencyStates: states.consistency,
            descriptions: states.descriptions,
            hasNotSavedEdits: hasNotSavedEdits.value,
            cb,
          });
        },
        getFuncCallCustomExports: (fc: DG.FuncCall) => {
          return Utils.getCustomExports(fc.func).map(x => x.name);
        },
        runFuncCallCustomExport: async (fc: DG.FuncCall, uuid: string, exportName: string) => {
          const exports =  Utils.getCustomExports(fc.func);
          const item = exports.find(x => x.name === exportName);
          if (!item)
            throw new Error(`No export named ${exportName} is defined for ${fc.func.nqName}`);
          const res = await DG.Func.byName(fc.func.nqName).apply({
            startDownload: false,
            funcCall: fc,
            validationState: states.validations?.[uuid],
            consistencyState: states.consistency?.[uuid],
            isOutputOutdated: states.calls?.[uuid]?.isOutputOutdated,
            runError: states.calls?.[uuid]?.runError,
          });
          return res;
        },
      };
      await exportData.handler(treeState.value!, utils);
    }

    ////
    // chosen step state
    ////

    const chosenStepState = Vue.computed(() => chosenStep.value?.state);

    const isRunDisabled = Vue.computed(() => {
      if (!chosenStepUuid.value)
        return true;
      const callState = states.calls[chosenStepUuid.value];
      if (!callState)
        return true;
      return (!callState.isRunnable || callState.isRunning || treeMutationsLocked.value || chosenStepState.value?.isReadonly);
    });

    const isOutputOutdated = Vue.computed(() => {
      if (!chosenStepUuid.value)
        return false;
      const callState = states.calls[chosenStepUuid.value];
      if (!callState)
        return false;
      return callState.isOutputOutdated;
    });

    Vue.watch(chosenStepState, (state) => {
      if (!state)
        return;
      if (isFuncCallState(state)) rfvHidden.value = false;
      if (!isFuncCallState(state)) pipelineViewHidden.value = false;
    });

    const isRootChoosen = Vue.computed(() => {
      return (!!chosenStepState.value?.uuid) && chosenStepState.value?.uuid === treeState.value?.uuid;
    });

    const menuActions = Vue.computed(() => {
      if (!treeState.value || !chosenStepUuid.value)
        return {};
      const globalActions = getRelevantGlobalActions(treeState.value, chosenStepUuid.value);
      const currentStepActions = chosenStepState.value?.actions?.filter(action => action.position === 'menu') ?? [];
      const actions = [...globalActions, ...currentStepActions];
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
        if (action.position === 'buttons')
          acc.push(action);

        return acc;
      }, [] as ViewAction[]);
    });

    ////
    // Navigation tree related
    ////

    let oldClosed = new Set<string>();

    Vue.watch(treeState, (_nextState, oldState) => {
      if (oldState && treeInstance.value) {
        const oldStats = treeInstance.value.statsFlat as AugmentedStat[];
        oldClosed = oldStats.reduce((acc, stat) => {
          if (!stat.open)
            acc.add(stat.data.uuid);
          return acc;
        }, new Set<string>());
      }
    }, {immediate: true});

    const restoreOpenedNodes = (stat: AugmentedStat) => {
      if (oldClosed.has(stat.data.uuid))
        stat.open = false;
      return stat;
    };

    const isTreeReportable = Vue.computed(() => {
      return Object.values(states.calls)
        .map((state) => state?.isOutputOutdated)
        .every((isOutdated) => isOutdated === false);
    });

    const isEachDraggable = (stat: AugmentedStat) => {
      return stat.parent && !stat.parent.data.isReadonly &&
        (isParallelPipelineState(stat.parent.data) || isSequentialPipelineState(stat.parent.data)) &&
        !stat.parent.data.stepTypes.find((item) => item.configId === stat.data.configId && item.disableUIDragging);
    };

    const isEachDroppable = (stat: AugmentedStat) => {
      const draggedStep = dragContext?.startInfo?.dragNode as AugmentedStat | undefined;
      return (isParallelPipelineState(stat.data) || isSequentialPipelineState(stat.data)) &&
        (!draggedStep || stat.data.uuid === draggedStep.parent?.data.uuid);
    };

    const onAfterDrop = () => {
      const draggedStep = dragContext.startInfo?.dragNode as AugmentedStat | undefined;
      if (draggedStep) {
        const oldIndex = dragContext.startInfo.indexBeforeDrop;
        const newIndex = dragContext.startInfo.indexBeforeDrop < dragContext.targetInfo.indexBeforeDrop ?
          dragContext.targetInfo.indexBeforeDrop - 1:
          dragContext.targetInfo.indexBeforeDrop;
        if (oldIndex !== newIndex)
          moveStep(draggedStep.data.uuid, newIndex);
      }
    };

    const isDeletable = (stat: AugmentedStat) => {
      return !!stat.parent && !stat.parent.data.isReadonly &&
        (isParallelPipelineState(stat.parent.data) || isSequentialPipelineState(stat.parent.data)) &&
        !stat.parent.data.stepTypes.find((item) => item.configId === stat.data.configId && item.disableUIRemoving);
    };

    ////
    // additional
    ////

    const isUserDeveloper = Vue.ref(false);
    Vue.onMounted(async () => {
      // Workaround till JS API is not ready: https://reddata.atlassian.net/browse/GROK-14159
      const userGroups = (await(await fetch(`${window.location.origin}/api/groups/all_parents`)).json() as DG.Group[]);

      if (userGroups.find((group) => group.friendlyName === DEVELOPERS_GROUP))
        isUserDeveloper.value = true;
    });

    ////
    // render
    ////

    return () => (
      Vue.withDirectives(<div class='w-full h-full'>
        <RibbonPanel view={currentView.value}>
          <IconFA
            name='folder-tree'
            tooltip={treeHidden.value ? 'Show tree': 'Hide tree'}
            onClick={() => treeHidden.value = !treeHidden.value }
          />
          { isUserDeveloper.value && <IconFA
            name='bug'
            tooltip={inspectorHidden.value ? 'Show inspector': 'Hide inspector'}
            onClick={() => inspectorHidden.value = !inspectorHidden.value }
          /> }
          {isTreeReady.value &&
            treeState.value &&
            (hasSubtreeFixableInconsistencies(treeState.value, states.calls, states.consistency) ?
              <IconFA
                name='sync'
                tooltip={'Update tree with consistent values'}
                style={{'padding-right': '3px'}}
                onClick={() => runSubtreeWithConfirm(treeState.value!.uuid, true)}
              />:
              <IconFA
                name='forward'
                tooltip={'Run ready steps'}
                style={{'padding-right': '3px'}}
                onClick={() => runSequence(treeState.value!.uuid, false)}
              />
            )}
          {isTreeReady.value && <IconFA
            name='save'
            tooltip={'Save current state of model'}
            style={{'padding-right': '3px'}}
            onClick={saveEntireModelState}
          /> }
          {isTreeReady.value && showReturn.value && <IconFA
            name='check'
            tooltip={'Confim data'}
            style={{'padding-right': '3px'}}
            onClick={onReturnClicked}
          />
          }
        </RibbonPanel>
        {isTreeReady.value && isTreeReportable.value &&
          <RibbonMenu groupName='Export' view={currentView.value}>
            {
              exports.value.map((exportData) =>
                <span onClick={() => exportHandler(exportData)}>
                  <div> {exportData.friendlyName ?? exportData.id} </div>
                </span>,
              )
            }
          </RibbonMenu>
        }
        {(reportBugUrl.value || reqFeatureUrl.value) &&
          <RibbonMenu groupName='Feedback' view={currentView.value}>
            {
              reportBugUrl.value &&
              <span onClick={() => window.open(reportBugUrl.value, '_blank')}>
                <div>Report a bug</div>
              </span>
            }
            {
              reqFeatureUrl.value &&
              <span onClick={() => window.open(reqFeatureUrl.value, '_blank')}>
                <div>Request a feature</div>
              </span>
            }
          </RibbonMenu>
        }
        { menuActions.value && Object.entries(menuActions.value).map(([category, actions]) =>
          <RibbonMenu groupName={category} view={currentView.value}>
            {
              actions.map((action) => Vue.withDirectives(<span onClick={() => runActionWithConfirmation(action.uuid)}>
                <div> { action.icon && <IconFA name={action.icon} style={{width: '15px', display: 'inline-block', textAlign: 'center'}}/> } { action.friendlyName ?? action.id } </div>
              </span>, [[tooltip, action.description]]))
            }
          </RibbonMenu>)
        }
        <DockManager class='block h-full'
          style={{overflow: 'hidden !important'}}
          key={providerFunc.value?.id}
          onPanelClosed={handlePanelClose}
        >
          { !inspectorHidden.value &&
            <Inspector
              key="inspector"
              treeState={treeState.value}
              config={config.value}
              logs={logs.value}
              links={links.value}
              ref={inspectorInstance}
              dock-spawn-title='Inspector'
              class='h-full overflow-scroll'
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio='0.2'
            ></Inspector>
          }
          {
            treeState.value && !treeHidden.value ?
              Vue.withDirectives(<Draggable
                class="ui-div mtl-tree p-2 overflow-scroll h-full"
                style={{paddingLeft: '25px'}}
                key="navigation"

                dock-spawn-title='Steps'
                dock-spawn-panel-icon='folder-tree'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-z-index={2}

                rootDroppable={false}
                treeLine
                childrenKey='steps'
                nodeKey={(stat: AugmentedStat) => stat.data.uuid}
                statHandler={restoreOpenedNodes}

                ref={treeInstance}
                modelValue={[treeState.value]}

                eachDraggable={isEachDraggable}
                eachDroppable={isEachDroppable}
                onAfter-drop={onAfterDrop}
                onClick:node={(stat) => {
                  chosenStepUuid.value = stat.data.uuid;

                  if (isFuncCallState(stat.data)) rfvHidden.value = false;
                  if (!isFuncCallState(stat.data)) pipelineViewHidden.value = false;
                }}
                onClose:node={(stat) => oldClosed.add(stat.data.uuid)}
                onOpen:node={(stat) => oldClosed.delete(stat.data.uuid)}
              >
                {
                  ({stat}: {stat: AugmentedStat}) =>
                    (
                      <TreeNode
                        stat={stat}
                        callState={states.calls[stat.data.uuid]}
                        validationStates={states.validations[stat.data.uuid]}
                        consistencyStates={states.consistency[stat.data.uuid]}
                        descriptions={states.descriptions[stat.data.uuid]}
                        style={{
                          'background-color': stat.data.uuid === chosenStepUuid.value ? '#f2f2f5' : null,
                        }}
                        isDraggable={treeInstance.value?.isDraggable(stat)}
                        isDroppable={treeInstance.value?.isDroppable(stat)}
                        isDeletable={isDeletable(stat)}
                        isReadonly={stat.data.isReadonly}
                        hasInconsistentSubsteps={!!hasSubtreeFixableInconsistencies(stat.data, states.calls, states.consistency)}
                        onAddNode={({itemId, position}) => addStep(stat.data.uuid, itemId, position)}
                        onRemoveNode={() => removeStep(stat.data.uuid)}
                        onToggleNode={() => stat.open = !stat.open}
                        onRunSubtree={(startUuid, rerunWithConsistent) => runSubtreeWithConfirm(startUuid, rerunWithConsistent)}
                        onRunStep={(uuid) => runStep(uuid)}
                        onSaveStep={(uuid) => saveSubTreeState(uuid)}
                      />
                    )
                }
              </Draggable>, [[ifOverlapping, treeMutationsLocked.value, 'Locked...']]): null
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
                skipInit={true}
                onUpdate:funcCall={onFuncCallChange}
                onActionRequested={runActionWithConfirmation}
                onConsistencyReset={(ioName) => consistencyReset(chosenStepUuid.value!, ioName)}
                dock-spawn-title='Step review'
                ref={rfvRef}
                view={currentView.value}
              >
                {{
                  navigation: ({runLabel, allowRerun}: {runLabel: string, allowRerun: boolean}) => (
                    <>
                      {
                        <Button onClick={goBack} style={{'margin-right': 'auto'}}>
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
                      {
                        (isOutputOutdated.value || allowRerun) &&
                          <BigButton
                            isDisabled={isRunDisabled.value}
                            onClick={() => runStep(chosenStepUuid.value!)}
                          >
                            {isOutputOutdated.value ? runLabel: 'Rerun'}
                          </BigButton>
                      }
                      {
                        !isOutputOutdated.value && !hasInconsistencies(states.consistency[chosenStepUuid.value!]) &&
                          <BigButton
                            onClick={goNextStep}
                          >
                             Next
                          </BigButton>
                      }
                      { hasInconsistencies(states.consistency[chosenStepUuid.value!]) &&
                        <BigButton
                          onClick={() => runSequence(chosenStepUuid.value!, true)}
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
            !pipelineViewHidden.value && chosenStepUuid.value && chosenStepState.value && !isFuncCallState(chosenStepState.value) &&
            <PipelineView
              funcCall={chosenStepState.value.nqName ?
                DG.Func.byName(chosenStepState.value.nqName!).prepare() :
                undefined
              }
              key={chosenStepUuid.value!}
              state={chosenStepState.value}
              uuid={chosenStepUuid.value}
              isRoot={isRootChoosen.value}
              buttonActions={buttonActions.value}
              onActionRequested={runActionWithConfirmation}
              dock-spawn-title='Step sequence review'
              onProceedClicked={onPipelineProceed}
              onBackClicked={goBack}
              onUpdate:funcCall={onPipelineFuncCallUpdate}
              onAddNode={({itemId, position}) => {
                if (chosenStepUuid.value)
                  addStep(chosenStepUuid.value, itemId, position);
              }}
              ref={pipelineViewRef}
              view={currentView.value}
            />
          }
        </DockManager>
      </div>, [[ifOverlapping, isGlobalLocked.value || !treeState.value]])
    );
  },
});
