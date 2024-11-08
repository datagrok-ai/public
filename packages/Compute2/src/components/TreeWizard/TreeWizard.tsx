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
import {BehaviorSubject, of} from 'rxjs';
import {PipelineView} from '../PipelineView/PipelineView';
import {useUrlSearchParams} from '@vueuse/core';
import {Inspector} from '../Inspector/Inspector';
import {switchMap} from 'rxjs/operators';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {findTreeNode, findTreeNodeParrent, makeMergedItems, reportStep} from '../../utils';


export const TreeWizard = Vue.defineComponent({
  name: 'TreeWizard',
  props: {
    providerFunc: {type: String, required: true},
  },
  setup(props) {
    const driver = new Driver();
    const treeState = Vue.shallowRef<PipelineState | undefined>(undefined);
    const logs = useObservable(driver.logger.logs$);
    const config = useObservable(driver.currentConfig$);
    const links = useObservable(driver.currentLinks$);

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

    useSubscription((driver.currentState$).subscribe((s) => {
      if (treeInstance.value) {
        const oldStats = treeInstance.value!.statsFlat as AugmentedStat[];
        oldClosed = oldStats.reduce((acc, stat) => {
          if (!stat.open)
            acc.push(stat.data.uuid);
          return acc;
        }, [] as string[]);
      }
      treeState.value = s;
    }));

    const treeMutationsLocked = useSubject(driver.treeMutationsLocked$);
    const isGlobalLocked = useSubject(driver.globalROLocked$);

    const states = Vue.reactive({
      calls: {} as Record<string, FuncCallStateInfo | undefined>,
      validations: {} as Record<string, Record<string, ValidationResult> | undefined>,
      consistency: {} as Record<string, Record<string, ConsistencyInfo> | undefined>,
      meta: {} as Record<string, Record<string, BehaviorSubject<any>> | undefined>
    });

    useSubscription(driver.currentCallsState$.pipe(
      switchMap((data) => {
        states.calls = {};
        return makeMergedItems(data);
      }),
    ).subscribe(([k, val]) => {
      states.calls[k] = Object.freeze(val);
    }));

    useSubscription(driver.currentValidations$.pipe(
      switchMap((data) => {
        states.validations = {};
        return makeMergedItems(data);
      }),
    ).subscribe(([k, val]) => {
      states.validations[k] = Object.freeze(val);
    }));

    useSubscription(driver.currentConsistency$.pipe(
      switchMap((data) => {
        states.consistency = {};
        return makeMergedItems(data);
      }),
    ).subscribe(([k, val]) => {
      states.consistency[k] = Object.freeze(val);
    }));

    useSubscription(driver.currentMeta$.pipe(
      switchMap((data) => {
        states.meta = {};
        return makeMergedItems(data);
      }),
    ).subscribe(([k, val]) => {
      states.meta[k] = Object.freeze(val);
    }));

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

    Vue.onMounted(() => {
      initPipeline(props.providerFunc);
    });

    Vue.onUnmounted(() => {
      driver.close();
    });

    const initPipeline = (provider: string) => {
      driver.sendCommand({event: 'initPipeline', provider});
    };

    const loadPipeline = (funcCallId: string) => {
      driver.sendCommand({event: 'loadPipeline', funcCallId});
    };

    const loadAndReplaceNestedPipeline = (parentUuid: string, dbId: string, itemId: string, position: number) => {
      driver.sendCommand({event: 'loadDynamicItem', parentUuid, dbId, itemId, position, readonly: true, isReplace: true});
    }

    const savePipeline = () => {
      driver.sendCommand({event: 'savePipeline'})
    };

    const runStep = async (uuid: string) => {
      driver.sendCommand({event: 'runStep', uuid});
    };

    const addStep = (parentUuid: string, itemId: string, position: number) => {
      driver.sendCommand({event: 'addDynamicItem', parentUuid, itemId, position});
    };

    const removeStep = (uuid: string) => {
      driver.sendCommand({event: 'removeDynamicItem', uuid});
    };

    const moveStep = (uuid: string, position: number) => {
      driver.sendCommand({event: 'moveDynamicItem', uuid, position});
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
        const newIndex = dragContext.startInfo.indexBeforeDrop < dragContext.targetInfo.indexBeforeDrop ?
          dragContext.targetInfo.indexBeforeDrop - 1:
          dragContext.targetInfo.indexBeforeDrop;
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
                dock-spawn-z-index={102}

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
                        onClick={() => {
                          chosenStepUuid.value = stat.data.uuid;
                        }}
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
                dock-spawn-title='Step review'
                ref={rfvRef}
              />
          }
          {
            chosenStepState.value &&
            !isFuncCallState(chosenStepState.value) && chosenStepState.value.provider &&
            <PipelineView
              funcCall={DG.Func.byName(chosenStepState.value.nqName!).prepare()}
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
            />
          }
        </DockManager> }
      </div>, [[ifOverlapping, isGlobalLocked.value]])
    );
  },
});
