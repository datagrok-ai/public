import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {DockManager, IconFA, ifOverlapping, RibbonMenu, RibbonPanel} from '@datagrok-libraries/webcomponents-vue';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {useObservable, useSubject, useSubscription, useExtractedObservable} from '@vueuse/rxjs';
import {
  isFuncCallState, isParallelPipelineState,
  isSequentialPipelineState, PipelineState,
  StepFunCallState,
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
import {findTreeNode, findTreeNodeParrent, reportStep} from '../../utils';
import {useReactiveTreeDriver} from '../../composables/use-reactive-tree-driver';

export const TreeWizard = Vue.defineComponent({
  name: 'TreeWizard',
  props: {
    providerFunc: {type: String, required: true},
  },
  setup(props) {

    // TODO: handle providerFunc changes, not necessary as of now
    const {
      treeMutationsLocked,
      isGlobalLocked,
      treeState,
      states,
      logs,
      config,
      links,
      //
      loadPipeline,
      loadAndReplaceNestedPipeline,
      savePipeline,
      runStep,
      runAction,
      addStep,
      removeStep,
      moveStep,
    } = useReactiveTreeDriver(props.providerFunc);


    const chosenStepUuid = Vue.ref<string | undefined>(undefined);

    const searchParams = useUrlSearchParams('history');

    const handleActivePanelChanged = async (newPanel: string | null, prevPanel: string | null) => {
      if (prevPanel === 'Steps') return;

      if (prevPanel === 'Step review') {
        rfvRef.value?.savePersonalState();
      }

      if (newPanel === 'Step review') {
        await Vue.nextTick();
        setTimeout(() => {
          rfvRef.value?.loadPersonalLayout();
        }, 50)
      }
    };

    Vue.watch(treeState, () => {
      if (!treeState.value) return;

      // Getting inital URL user entered with
      const startUrl = new URL(grok.shell.startUri);
      const stepId = startUrl.searchParams.get('stepId');

      if (!stepId) return;

      chosenStepUuid.value = stepId;
    });

    const chosenStepState = Vue.computed(() => {
      if (!chosenStepUuid.value || !treeState.value) return null;

      return findTreeNode(chosenStepUuid.value, treeState.value);
    });

    const isRootChoosen = Vue.computed(() => {
      return (!!chosenStepState.value?.uuid) && chosenStepState.value?.uuid === treeState.value?.uuid;
    })

    Vue.watch(chosenStepUuid, (newStepId) => {
      if (newStepId)
        searchParams['stepId'] = newStepId;
    });

    const treeInstance = Vue.ref(null as InstanceType<typeof Draggable> | null);

    let oldClosed = [] as string[];

    Vue.watch(treeState, (_nextState, oldState) => {
      if (oldState && treeInstance.value) {
        const oldStats = treeInstance.value.statsFlat as AugmentedStat[];
        oldClosed = oldStats.reduce((acc, stat) => {
          if (!stat.open)
            acc.push(stat.data.uuid);
          return acc;
        }, [] as string[]);
      }
    }, { immediate: true });

    const isTreeReportable = Vue.computed(() => {
      return Object.values(states.calls)
        .map((state) => state?.isOutputOutdated)
        .every((isOutdated) => isOutdated === false);
    });

    const restoreOpenedNodes = (stat: AugmentedStat) => {
      if (oldClosed.includes(stat.data.uuid))
        stat.open = false;
      return stat;
    };

    const treeHidden = Vue.ref(false);
    const inspectorHidden = Vue.ref(true);
    const inspectorInstance = Vue.ref(null as InstanceType<typeof Inspector> | null);
    const rfvRef = Vue.ref(null as InstanceType<typeof RichFunctionView> | null);

    const handlePanelClose = (el: HTMLElement) => {
      if (el === treeInstance.value?.$el.parentElement) treeHidden.value = true;
      if (el === inspectorInstance.value?.$el) inspectorHidden.value = true;
    };

    const isEachDraggable = (stat: AugmentedStat) => {
      return (stat.parent && !stat.parent.data.isReadonly &&
        (isParallelPipelineState(stat.parent.data) || isSequentialPipelineState(stat.parent.data))
      ) ?? false
    }

    const isEachDroppable = (stat: AugmentedStat) => {
      const draggedStep = dragContext?.startInfo?.dragNode as AugmentedStat | undefined;
      return (isParallelPipelineState(stat.data) || isSequentialPipelineState(stat.data)) &&
        (!draggedStep || stat.data.uuid === draggedStep.parent?.data.uuid);
    }

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
    }

    const isDeletable = (stat: AugmentedStat) => {
      return !!stat.parent && !stat.parent.data.isReadonly &&
        (isParallelPipelineState(stat.parent.data) || isSequentialPipelineState(stat.parent.data))
    }

    return () => (
      Vue.withDirectives(<div class='w-full h-full'>
        <RibbonPanel>
          <IconFA
            name='folder-tree'
            tooltip={treeHidden.value ? 'Show tree': 'Hide tree'}
            onClick={() => treeHidden.value = !treeHidden.value }
          />
          <IconFA
            name='bug'
            tooltip={inspectorHidden.value ? 'Show inspector': 'Hide inspector'}
            onClick={() => inspectorHidden.value = !inspectorHidden.value }
          />

          {treeState.value && isTreeReportable.value && <IconFA
            name='arrow-to-bottom'
            tooltip='Report all steps'
            onClick={async () => reportStep(treeState.value) }
          /> }
        </RibbonPanel>
        <RibbonMenu groupName='State'>
          <span onClick={savePipeline}>
            <IconFA name='save' style={{'padding-right': '3px'}}/>
            <span> Save </span>
          </span>
        </RibbonMenu>
        {treeState.value &&
        <DockManager class='block h-full'
          onPanelClosed={handlePanelClose}
          onUpdate:activePanelTitle={handleActivePanelChanged}
        >
          { !inspectorHidden.value &&
            <Inspector
              treeState={treeState.value}
              config={config.value}
              logs={logs.value}
              links={links.value}
              ref={inspectorInstance}
              dock-spawn-title='Inspector'
              class='h-full overflow-scroll'
            ></Inspector>
          }
          {
            treeState.value && !treeHidden.value ?
              Vue.withDirectives(<Draggable
                class="ui-div mtl-tree p-2 overflow-scroll h-full"
                style={{paddingLeft: '25px'}}

                dock-spawn-title='Steps'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.3}
                dock-spawn-z-index={51}

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
                  chosenStepUuid.value = stat.data.uuid
                }}
              >
                {
                  ({stat}: {stat: AugmentedStat}) =>
                    (
                      <TreeNode
                        stat={stat}
                        callState={states.calls[stat.data.uuid]}
                        validationStates={states.validations[stat.data.uuid]}
                        consistencyStates={states.consistency[stat.data.uuid]}
                        style={{
                          'background-color': stat.data.uuid === chosenStepUuid.value ? '#f2f2f5' : null,
                        }}
                        isDraggable={treeInstance.value?.isDraggable(stat)}
                        isDroppable={treeInstance.value?.isDroppable(stat)}
                        isDeletable={isDeletable(stat)}
                        isReadonly={stat.data.isReadonly}
                        onAddNode={({itemId, position}) => {
                          addStep(stat.data.uuid, itemId, position);
                        }}
                        onRemoveNode={() => removeStep(stat.data.uuid)}
                        onToggleNode={() => stat.open = !stat.open}
                      />
                    )
                }
              </Draggable>, [[ifOverlapping, treeMutationsLocked.value, 'Locked...']]): null
          }
          {
            chosenStepState.value &&
            isFuncCallState(chosenStepState.value) && chosenStepState.value.funcCall &&
              <RichFunctionView
                class='overflow-hidden'
                funcCall={chosenStepState.value.funcCall!}
                callState={chosenStepUuid.value ? states.calls[chosenStepUuid.value] : undefined}
                callMeta={chosenStepUuid.value ? states.meta[chosenStepUuid.value] : undefined}
                viewersHook={chosenStepState.value.viewersHook}
                validationStates={states.validations[chosenStepState.value.uuid]}
                isReadonly={chosenStepState.value.isReadonly}
                isTreeLocked={treeMutationsLocked.value}
                onUpdate:funcCall={(call) => (chosenStepState.value as StepFunCallState).funcCall = call}
                onRunClicked={() => runStep(chosenStepState.value!.uuid)}
                onActionRequested={runAction}
                dock-spawn-title='Step review'
                ref={rfvRef}
              />
          }
          {
            chosenStepState.value &&
            !isFuncCallState(chosenStepState.value) && chosenStepState.value.provider &&
            <PipelineView
              funcCall={DG.Func.byName(chosenStepState.value.nqName!).prepare()}
              state={chosenStepState.value}
              isRoot={isRootChoosen.value}
              dock-spawn-title='Step sequence review'
              onProceedClicked={() => {
                if (chosenStepState.value && !isFuncCallState(chosenStepState.value))
                  chosenStepUuid.value = chosenStepState.value.steps[0].uuid;
              }}
              onUpdate:funcCall={
                (newCall) => {
                  if (isRootChoosen.value) {
                    loadPipeline(newCall.id);
                  } else if (chosenStepState.value && chosenStepUuid.value && treeState.value && !isFuncCallState(chosenStepState.value)) {
                    const parent = findTreeNodeParrent(chosenStepUuid.value, treeState.value);
                    if (parent && !isFuncCallState(parent)) {
                      const position = parent.steps.findIndex(step => step.uuid === chosenStepUuid.value);
                      loadAndReplaceNestedPipeline(parent.uuid, newCall.id, chosenStepState.value.configId, position);
                    }
                  }
                }
              }
              onAddNode={({itemId, position}) => {
                addStep(chosenStepState.value!.uuid, itemId, position);
              }}
            />
          }
        </DockManager> }
      </div>, [[ifOverlapping, isGlobalLocked.value]])
    );
  },
});
