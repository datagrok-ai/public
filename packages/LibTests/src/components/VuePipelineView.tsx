import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {defineComponent, PropType, ref, shallowRef, triggerRef} from 'vue';
import {SplitH} from '@datagrok-libraries/webcomponents-vue/src';
import {Draggable, OpenIcon} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import {VueRichFunctionView} from './VueRichFunctionView';

type TreeNode = {
  text: string,
  children: TreeNode[]
}

type State = {
  checked: boolean,
  open: boolean,
  children: any[],
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

    return () => (
      <SplitH resize={true} style={{height: '100%', display: 'block', padding: '8px'}}>
        <div>
          <h2> Model name </h2>
          <Draggable v-model={props.treeData} treeLine> 
            { 
              ({stat, node}: {stat: State, node: TreeNode}) => {
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
                  console.log(currentFuncCall.value.inputs['initTemp']);
                  triggerRef(currentFuncCall);
                };
                return (
                  <div style={{display: 'flex'}} onClick={onNodeClick}>
                    { stat.children.length ? openIcon : null }
                    <input
                      class="mtl-checkbox mtl-mr"
                      type="checkbox"
                      v-model={stat.checked}
                      title='check'
                    />
                    <span class="mtl-ml">{node.text}</span>
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
