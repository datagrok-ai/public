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
  funcCall: DG.FuncCall | string,
  text: string
}

type AugmentedStat = Stat<Data> & {
  isHovered: boolean
  status: Status
};

const statusToIcon = {
  ['locked']: 'lock',
  [`didn't run`]: 'circle',
  ['running']: 'hourglass-half',
  ['succeeded']: 'check-circle',
  ['partially succeeded']: 'dot-circle',
  ['failed']: 'times-circle',
} as Record<Status, string>;

const statusToColor = {
  ['locked']: 'black',
  [`didn't run`]: 'yellow',
  ['running']: 'blue',
  ['succeeded']: 'green',
  ['partially succeeded']: 'green',
  ['failed']: 'red',
} as Record<Status, string>;

const statusToTooltip = {
  ['locked']: 'This step is currently locked',
  [`didn't run`]: `This step didn't run yet`,
  ['running']: 'Running...',
  ['succeeded']: 'Run succeeded',
  ['partially succeeded']: 'Partially succeeded',
  ['failed']: 'Run failed',
} as Record<Status, string>;

type Status = 'locked' | `didn't run` | 'running' | 'succeeded' | 'failed' | 'partially succeeded';

const getCall = (funcCall: DG.FuncCall | string, params?: {
  a?: number,
  b?: number,
}) => {
  const defaultParams = {
    a: 2,
    b: 3,
  };

  return funcCall instanceof DG.FuncCall ?
    funcCall.func.prepare({
      ...defaultParams,
      ...params,
    }) : DG.Func.byName(funcCall).prepare({
      ...defaultParams,
      ...params,
    });
};

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
            style={{paddingLeft: '20px'}}
            treeLine
          > 
            { 
              ({stat, node}: {stat: AugmentedStat, node: TreeNode}) => {
                const openIcon = <OpenIcon
                  open={stat.open}
                  class="mtl-mr"
                  //@ts-ignore
                  onClick={(e) => {stat.open = !stat.open; e.stopPropagation();}}
                />;

                const progressIcon = (status: Status) =>{ 
                  return <IconFA 
                    name={statusToIcon[status]}
                    animation={status === `running` ? 'spin': null}
                    tooltip={statusToTooltip[status] ?? null}
                    style={{
                      color: statusToColor[status],
                      alignSelf: 'center',
                      left: '-16px',
                      position: 'absolute',
                    }} 
                  />;
                };

                const onNodeClick = async () => {
                  stat.status = 'running';
                  if (stat.data.funcCall) {
                    try {
                      const call = await getCall(stat.data.funcCall).call();
                      currentFuncCall.value = call;
                      stat.status = 'succeeded';
                    } catch {
                      stat.status = 'failed';
                    }
                  
                    triggerRef(currentFuncCall);
                  } else {
                    await Promise.all(stat.children.map(async (child) => {
                      try {
                        child.status = 'running';
                        const call = await getCall(child.data.funcCall).call();
                        child.data.funcCall = call;
                        child.status = 'succeeded';

                        return Promise.resolve(child.status);
                      } catch (e) {
                        child.status = 'failed';
                        return Promise.reject(e);
                      }
                    })).then(() => {
                      stat.status = 'succeeded';
                    }).catch(() => {
                      stat.status = 'failed';
                    });                      
                  }
                };

                return (
                  <div style={{display: 'flex', padding: '4px 0px', width: '100%'}}
                    onClick={onNodeClick}
                    onMouseover={() => stat.isHovered = true} 
                    onMouseleave={() => stat.isHovered = false} 
                    onDragstart={() => stat.isHovered = false}
                  >
                    { progressIcon(stat.status ?? `didn't run`) }
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
                      onClick={(e: Event) => {tree.value!.remove(stat); e.stopPropagation();}}
                    />: null }
                    { tree.value?.isDroppable(stat) && stat.isHovered ? <IconFA 
                      name='plus' 
                      style={{paddingLeft: '4px'}}
                      onClick={(e) => {
                        tree.value!.add({
                          text: `${node.text.includes('phase') ? 'Phase': 'Day'} ${(stat.children.length + 1).toString()}`,
                        }, 
                        stat, stat.children.length);
                        e.stopPropagation();
                      }}
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
