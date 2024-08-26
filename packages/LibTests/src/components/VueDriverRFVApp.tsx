import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BigButton, Button, InputForm, SplitH} from '@datagrok-libraries/webcomponents-vue/src';
import {defineComponent, KeepAlive, onUnmounted, ref, shallowRef, triggerRef, watch} from 'vue';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {useSubscription} from '@vueuse/rxjs';
import {isFuncCallState, isParallelPipelineState, isSequentialPipelineState, PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {RichFunctionView} from './RFV/RichFunctionView';
import {TreeNode} from './TreeWizard/TreeNode';
import {Draggable} from '@he-tree/vue';
import {AugmentedStat, HueTree} from './TreeWizard/types';
import {dragContext} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';

export const VueDriverRFVApp = defineComponent({
  name: 'VueDriverRFVApp',
  setup() {
    const driver = new Driver();
    const isLocked = ref(false);
    const treeState = shallowRef<PipelineState | undefined>(undefined);

    const currentFuncCall = shallowRef<DG.FuncCall | undefined>(undefined);

    const treeInstance = ref(null as HueTree | null);

    const changeFunccall = (newCall: DG.FuncCall) => {
      isVisibleRfv.value = true;

      currentFuncCall.value = newCall;
    };

    const isVisibleRfv = ref(true);

    let oldClosed = [] as string[];

    useSubscription((driver.currentState$).subscribe((s) => {
      if (treeInstance.value) {
        console.log(treeInstance.value);
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

      ui.setUpdateIndicator(treeInstance.value!.$el, newVal);
    });

    return () => (
      <KeepAlive>
        <div style={{width: '100%', height: '100%'}}>
          <BigButton onClick={() => initPipeline('LibTests:MockProvider3')}>Init Pipeline</BigButton>
          
          { treeState.value ? <SplitH resize={true} style={{height: '100%', display: 'block', padding: '8px'}}>
            <Draggable 
              class="mtl-tree"
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
                      isDraggable={treeInstance.value?.isDraggable(stat)}
                      isDroppable={treeInstance.value?.isDroppable(stat)}
                      isDeletable={!!stat.parent && isParallelPipelineState(stat.parent.data)}
                      onAddNode={({itemId, position}) => {
                        addStep(stat.data.uuid, itemId, position);
                      }}
                      onRemoveNode={() => removeStep(stat.data.uuid)}
                      onClick={() => {
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
            </Draggable>
            <div>
              {
                isVisibleRfv.value && currentFuncCall.value && <RichFunctionView 
                  style={{height: '100%'}} 
                  funcCall={currentFuncCall.value}
                  onFuncCallChange={(chosenCall) => currentFuncCall.value=chosenCall}
                /> 
              }
            </div>
          </SplitH>: null }
        </div>
      </KeepAlive>
    );
  },
});
