import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {defineComponent, KeepAlive, PropType, ref, shallowRef, triggerRef} from 'vue';
import {SplitH} from '@datagrok-libraries/webcomponents-vue/src';
import {Draggable} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {RichFunctionView} from '../RFV/RichFunctionView';
import {Stat} from '@he-tree/vue/types/src/components/TreeProcessorVue';
import {TreeNode} from './TreeNode';
import {StepConfig, HueTree, AugmentedStat, Data} from './types';

const getInitialStatus = (stat: AugmentedStat) => {
  if (!stat.parent) return `didn't run`;
  const parent = stat.parent;

  if (parent.data.order === 'sequental') {
    if (parent.status === 'didn\'t run' && parent.children.at(0) === stat) 
      return `didn't run`;

    return 'locked';
  }

  if (parent.data.order === 'parallel' && parent.status !== `locked`) return `didn\'t run`;

  return `locked`;
};

export const TreeWizardView = defineComponent({
  props: {
    wrapperFunccall: {
      required: true,
      type: String,
    },
    treeData: {
      required: true,
      type: Object as PropType<StepConfig>,
    },
  },
  setup(props) {
    const initialName = 'Compute:ObjectCooling';
    const currentFuncCall = shallowRef(DG.Func.byName(initialName).prepare());

    const tree = ref(null as HueTree | null);

    const changeFunccall = (newCall: string) => {
      isVisibleRfv.value = true;

      currentFuncCall.value = DG.Func.byName(newCall).prepare();
      triggerRef(currentFuncCall);
    };

    const isVisibleRfv = ref(true);

    const freshStatHandler = (freshStat: AugmentedStat) => {
      freshStat.data.status = getInitialStatus(freshStat);
      
      return freshStat; 
    };

    return () => (
      <SplitH resize={true} style={{height: '100%', display: 'block', padding: '8px'}}>
        <Draggable 
          class="mtl-tree"
          rootDroppable={false}
          style={{paddingLeft: '20px'}}
          treeLine

          ref={tree} 
          v-model={props.treeData} 
          
          eachDraggable={(stat: Stat<Data>) => (stat.parent?.data.text.includes('phases') ?? false)}
          eachDroppable={(stat: Stat<Data>) => (stat.data.text.includes('Review'))}
          statHandler={freshStatHandler}
        > 
          { 
            ({stat, node}: {stat: AugmentedStat, node: StepConfig}) =>  
              (
                <TreeNode 
                  stat={stat}
                  node={node}
                  isDraggable={tree.value?.isDraggable(stat)}
                  isDroppable={tree.value?.isDroppable(stat)}
                  isDeletable={stat.parent && tree.value?.isDroppable(stat.parent)}
                  onAddNode={(text) => tree.value?.add({text}, stat, stat.children.length)}
                  onRemoveNode={() => tree.value?.remove(stat)}
                  onClick={() => {
                    if (stat.data.funcCall) 
                      changeFunccall(stat.data.funcCall); 
                    else 
                      isVisibleRfv.value = false;
                  }}
                />
              )
          }
        </Draggable>
        <div>
          {
            isVisibleRfv.value && <RichFunctionView 
              style={{height: '100%'}} 
              funcCall={currentFuncCall.value}
              onFuncCallChange={(chosenCall) => currentFuncCall.value=chosenCall}
            /> 
          }
        </div>
      </SplitH>
    );
  },
});
