import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BigButton, Button, DockManager, FoldableDialog, IconFA, InputForm, RibbonPanel, SplitH} from '@datagrok-libraries/webcomponents-vue';
import {defineComponent, KeepAlive, onUnmounted, ref, shallowRef, triggerRef, VNode, watch} from 'vue';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {useSubscription} from '@vueuse/rxjs';
import {isFuncCallState, isParallelPipelineState, isSequentialPipelineState, PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {RichFunctionView} from './RFV/RichFunctionView';
import {TreeNode} from './TreeWizard/TreeNode';
import {Draggable} from '@he-tree/vue';
import {AugmentedStat} from './TreeWizard/types';
import {dragContext} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {deepCopy} from '@datagrok-libraries/compute-utils/shared-utils/utils';

export const VueDriverRFVApp = defineComponent({
  name: 'VueDriverRFVApp',
  setup() {
    const driver = new Driver();
    const isLocked = ref(false);
    const treeState = shallowRef<PipelineState | undefined>(undefined);

    const currentFuncCall = shallowRef<DG.FuncCall | undefined>(undefined);

    const treeInstance = ref(null as InstanceType<typeof Draggable> | null);

    const changeFunccall = (newCall: DG.FuncCall) => {
      isVisibleRfv.value = true;

      currentFuncCall.value = newCall;
    };

    const isVisibleRfv = ref(true);

    let oldClosed = [] as string[];

    useSubscription((driver.currentState$).subscribe((s) => {
      if (treeInstance.value) {
        //@ts-ignore
        const oldStats = treeInstance.value!.statsFlat as AugmentedStat[];
        oldClosed = oldStats.reduce((acc, stat) => {
          if (!stat.open) 
            acc.push(stat.data.uuid); 
          return acc;
        }, [] as string[]);
      }
      treeState.value = s;
    }));

    const restoreOpenedNodes = (stat: AugmentedStat) => {
      if (oldClosed.includes(stat.data.uuid)) 
        stat.open = false;
      return stat;
    };

    useSubscription((driver.stateLocked$).subscribe((l) => isLocked.value = l));
    onUnmounted(() => {
      driver.close();
      console.log('VuewDriverTestApp driver closed');
    });

    const initPipeline = (provider: string) => {
      driver.sendCommand({event: 'initPipeline', provider});
    };

    const runStep = (uuid: string) => {
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

    watch(isLocked, (newVal) => {
      if (!treeInstance.value) return;

      ui.setUpdateIndicator(treeInstance.value.$el, newVal);
    });

    const chosenStepUuid = ref(null as string | null);

    const treeHidden = ref(false);
    const rfvRef = ref(null as InstanceType<typeof RichFunctionView> | null);

    const handlePanelClose = (el: HTMLElement) => {
      if (el === rfvRef.value?.$el) isVisibleRfv.value = false;
      if (el === treeInstance.value?.$el) treeHidden.value = true;
    };

    return () => (
      <div class='w-full h-full'>
        <RibbonPanel>
          <BigButton onClick={() => initPipeline('LibTests:MockProvider3')}>Init Pipeline</BigButton>
          <IconFA 
            name='folder-tree'
            tooltip={treeHidden.value ? 'Show tree': 'Hide tree'}
            onClick={() => treeHidden.value = !treeHidden.value } 
          />
        </RibbonPanel>
        <DockManager class='block h-full' onPanelClosed={handlePanelClose}>
          { treeState.value && !treeHidden.value ? <Draggable 
            class="ui-div mtl-tree p-2"
            {...{title: 'Steps'}}
            rootDroppable={false}
            treeLine
            childrenKey='steps'
            nodeKey={(stat: AugmentedStat) => stat.data.uuid}
            statHandler={restoreOpenedNodes}

            ref={treeInstance} 
            modelValue={[treeState.value]} 
        
            eachDraggable={(stat: AugmentedStat) =>
              (stat.parent && 
                (isParallelPipelineState(stat.parent.data) || isSequentialPipelineState(stat.parent.data))
              ) ?? false
            }
            eachDroppable={(stat: AugmentedStat) => 
              (isParallelPipelineState(stat.data) || isSequentialPipelineState(stat.data))}
            onAfter-drop={() => {
              const draggedStep = dragContext.startInfo.dragNode as AugmentedStat;
              moveStep(draggedStep.data.uuid, dragContext.targetInfo.indexBeforeDrop);
            }}
          > 
            { 
              ({stat}: {stat: AugmentedStat}) =>  
                (
                  <TreeNode 
                    stat={stat}
                    style={{'background-color': stat.data.uuid === chosenStepUuid.value ? '#f2f2f5' : null}}
                    isDraggable={treeInstance.value?.isDraggable(stat)}
                    isDroppable={treeInstance.value?.isDroppable(stat)}
                    isDeletable={!!stat.parent && isParallelPipelineState(stat.parent.data)}
                    onAddNode={({itemId, position}) => {
                      addStep(stat.data.uuid, itemId, position);
                    }}
                    onRemoveNode={() => removeStep(stat.data.uuid)}
                    onClick={() => {
                      chosenStepUuid.value = stat.data.uuid;
                      if (isFuncCallState(stat.data) && stat.data.funcCall) 
                        changeFunccall(stat.data.funcCall); 
                      else 
                        isVisibleRfv.value = false;
                    }}
                    onRunNode={() => runStep(stat.data.uuid)}
                    onToggleNode={() => stat.open = !stat.open}
                  />
                )
            }
          </Draggable>: null }
          
          {
            isVisibleRfv.value && currentFuncCall.value && <RichFunctionView 
              class='overflow-hidden'
              funcCall={currentFuncCall.value}
              onUpdate:funcCall={(chosenCall) => currentFuncCall.value = deepCopy(chosenCall)}
              {...{title: 'Step review'}}
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.8}
              ref={rfvRef}
            /> 
          }
        </DockManager>
      </div>
    );
  },
});
