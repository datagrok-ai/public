import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {defineComponent, PropType, ref, shallowRef, triggerRef, watch} from 'vue';
import {IconFA, SplitH} from '@datagrok-libraries/webcomponents-vue/src';
import {Draggable, OpenIcon} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {VueRichFunctionView} from './VueRichFunctionView';
import {Stat} from '@he-tree/vue/types/src/components/TreeProcessorVue';

type TreeNode = {
  text: string,
  children: TreeNode[]
}

type Data = {
  text: string
}

export const VuePipelineView = defineComponent({
  props: {
    wrapperFunccall: {
      required: true,
      type: Object as PropType<DG.FuncCall | string>,
    },
    treeData: {
      required: true,
      type: Object as PropType<TreeNode>,
    },
  },
  setup(props) {
    const currentFuncCall = shallowRef(props.wrapperFunccall);
    const tree = ref(null as {
      add: Function,
      remove: Function,
      isDraggable: Function,
      isDroppable: Function,
    } | null);

    return () => (
      <SplitH resize={true} style={{height: '100%', display: 'block', padding: '8px'}}>
        <div>
          <h2> Model name </h2>
          <Draggable 
            class="mtl-tree"
            ref={tree} 
            v-model={props.treeData} 
            eachDraggable={(stat: Stat<Data>) => (stat.data.text.includes('Phase'))}
            eachDroppable={(stat: Stat<Data>) => (stat.data.text.includes('Review'))}
            rootDroppable={false}
            treeLine
          > 
            { 
              ({stat, node}: {stat: Stat<Data>, node: TreeNode}) => {
                const openIcon = <OpenIcon
                  open={stat.open}
                  class="mtl-mr"
                  //@ts-ignore
                  onClick={() => {stat.open = !stat.open;}}
                />;

                const getCall = (params: {
                  ambTemp?: number,
                  initTemp?: number,
                  desiredTemp?: number,
                  area?: number, 
                  heatCap?: number,
                  heatTransferCoeff?: number,
                  simTime?: number,
                }) => {
                  const defaultParams = {
                    ambTemp: 22,
                    initTemp: 100,
                    desiredTemp: 30,
                    area: 0.06, 
                    heatCap: 4200,
                    heatTransferCoeff: 8.3,
                    simTime: 21600,
                  };

                  return props.wrapperFunccall instanceof DG.FuncCall ?
                    props.wrapperFunccall.func.prepare({
                      ...defaultParams,
                      ...params,
                    }) : DG.Func.byName(props.wrapperFunccall).prepare({
                      ...defaultParams,
                      ...params,
                    });
                };

                const onNodeClick = () => {
                  currentFuncCall.value = getCall({initTemp: Math.random()*70 + 30});
                  triggerRef(currentFuncCall);
                };
                return (
                  <div style={{display: 'flex', padding: '4px 0px', width: '100%'}}
                    onClick={onNodeClick}
                    onMouseover={() => stat.isHovered = true} 
                    onMouseleave={() => stat.isHovered = false} 
                    onDragstart={() => stat.isHovered = false}
                  >
                    { stat.children.length ? openIcon : null }
                    <span class="mtl-ml">{node.text}</span>
                    { tree.value?.isDraggable(stat) && stat.isHovered ? <IconFA 
                      name='grip-vertical' 
                      cursor='grab'
                      style={{paddingLeft: '4px'}}
                    />: null }
                    { stat.parent && tree.value?.isDroppable(stat.parent) && stat.isHovered ? <IconFA 
                      name='times' 
                      style={{paddingLeft: '4px'}}
                      onClick={() => tree.value!.remove(stat)}
                    />: null }
                    { tree.value?.isDroppable(stat) && stat.isHovered ? <IconFA 
                      name='plus' 
                      style={{paddingLeft: '4px'}}
                      onClick={() => tree.value!.add({
                        text: `${node.text.includes('phase') ? 'Phase': 'Day'} ${(stat.children.length + 1).toString()}`,
                      }, 
                      stat, stat.children.length)}
                    />: null }
                  </div>
                );
              }
            }
          </Draggable>
        </div>
        <VueRichFunctionView funcCall={currentFuncCall.value} style={{height: '100%'}}/>
      </SplitH>
    );
  },
});
