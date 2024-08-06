import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {defineComponent, PropType, ref, shallowRef} from 'vue';
import {SplitH} from '@datagrok-libraries/webcomponents-vue/src';
import {Draggable} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {RichFunctionView} from '../RFV/RichFunctionView';
import {Stat} from '@he-tree/vue/types/src/components/TreeProcessorVue';
import {TreeNode} from './TreeNode';
import {TreeNodeType, HueTree, AugmentedStat, Data} from './types';

export const TreeWizardView = defineComponent({
  props: {
    wrapperFunccall: {
      required: true,
      type: Object as PropType<DG.FuncCall | string>,
    },
    treeData: {
      required: true,
      type: Object as PropType<TreeNodeType>,
    },
  },
  setup(props) {
    const currentFuncCall = shallowRef(props.wrapperFunccall);
    const tree = ref(null as HueTree | null);

    const changeFunccall = (newCall: DG.FuncCall | string) => {
      isVisibleRfv.value = true;
      currentFuncCall.value = newCall;
    };

    const isVisibleRfv = ref(true);

    return () => (
      <SplitH resize={true} style={{height: '100%', display: 'block', padding: '8px'}}>
        <Draggable 
          class="mtl-tree"
          ref={tree} 
          v-model={props.treeData} 
          eachDraggable={(stat: Stat<Data>) => (stat.parent?.data.text.includes('phases') ?? false)}
          eachDroppable={(stat: Stat<Data>) => (stat.data.text.includes('Review'))}
          rootDroppable={false}
          style={{paddingLeft: '20px'}}
          treeLine
        > 
          { 
            ({stat, node}: {stat: AugmentedStat, node: TreeNodeType}) =>  
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
        <RichFunctionView style={{display: isVisibleRfv, height: '100%'}} funcCall={currentFuncCall.value}/> 
      </SplitH>
    );
  },
});
