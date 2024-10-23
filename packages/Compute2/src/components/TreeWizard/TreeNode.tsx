import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {AugmentedStat, Status} from './types';
import {ComboPopup, IconFA} from '@datagrok-libraries/webcomponents-vue';
import {OpenIcon} from '@he-tree/vue';
import {
  isFuncCallState, isParallelPipelineState, 
  isSequentialPipelineState, PipelineState, 
  PipelineStateParallel, PipelineStateSequential,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {useElementHover} from '@vueuse/core';
import {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';

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

const runSequentallly = async (children: AugmentedStat[]) => {
  for (const child of children) 
    await runByTree(child);
};

const runInParallel = async (children: AugmentedStat[]): Promise<void[]> => 
  Promise.all(children.map(async (child) => runByTree(child as AugmentedStat)));

const updateNextStepAndParent = (currentStat: AugmentedStat) => {
  const parent = currentStat.parent;
  if (!parent) return;

  const getNextSibling = (step?: AugmentedStat | null) => {
    const parent = step?.parent;
    if (!parent) return;

    return parent.children[
      parent.children.findIndex(
        (childStat) => childStat === step,
      ) + 1
    ];
  };

  const nextSibling = getNextSibling(currentStat);

  if (nextSibling) 
    nextSibling.status = `didn't run`;
  else {
    // parent.status = currentStat.data.status;
    const parentNextSibling = getNextSibling(parent);
    if (parentNextSibling)
      parentNextSibling.status = `didn't run`;
  }
};

const runByTree = async (currentStat: AugmentedStat) => {
  const nodeCall = isFuncCallState(currentStat.data) ? currentStat.data.funcCall: null;
  // currentStat.data.status = 'running';
  if (nodeCall) {
    try {
      await getCall(nodeCall).call();
      // currentStat.data.status = 'succeeded';
      
      updateNextStepAndParent(currentStat);

      // return Promise.resolve(currentStat.data.status);
      return Promise.resolve();
    } catch (e) {
      // currentStat.data.status = 'failed';

      return Promise.reject(e);
    }
  } else {  
    return (isParallelPipelineState(currentStat.data) ? 
      runInParallel(currentStat.children): 
      runSequentallly(currentStat.children)
    ).then(() => {
      // currentStat.data.status = 'succeeded';

      updateNextStepAndParent(currentStat);

      // return Promise.resolve(currentStat.data.status);
      return Promise.resolve();
    }).catch((e) => {
      // currentStat.data.status = 'failed';

      return Promise.reject(e);
    });                      
  }
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
  [`didn't run`]: 'blue',
  ['running']: 'blue',
  ['succeeded']: 'green',
  ['partially succeeded']: 'red',
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

const callStateToStatus = (callState: FuncCallStateInfo): Status => {
  if (callState.isRunning) return 'running';
  if (!callState.isOutputOutdated) return 'succeeded';
  
  return 'didn\'t run';
};

export const TreeNode = Vue.defineComponent({
  props: {
    stat: {
      type: Object as Vue.PropType<AugmentedStat>,
      required: true,
    },
    callState: {
      type: Object as Vue.PropType<FuncCallStateInfo>,
    },
    isDraggable: {
      type: Boolean,
    },
    isDroppable: {
      type: Boolean,
    },
    isDeletable: {
      type: Boolean,
    },
  },
  emits: {
    addNode: ({itemId, position}:{itemId: string, position: number}) => ({itemId, position}),
    removeNode: () => {},
    click: () => {},
    toggleNode: () => {},
  },
  setup(props, {emit}) {      
    const openIcon = () => <OpenIcon
      open={props.stat.open}
      class="mtl-mr"
      //@ts-ignore
      onClick={(e) => {emit('toggleNode'); e.stopPropagation();}}
    />;
    
    const progressIcon = (status: Status) =>{ 
      return <IconFA 
        name={statusToIcon[status]}
        animation={status === `running` ? 'spin': null}
        tooltip={statusToTooltip[status] ?? null}
        style={{
          color: statusToColor[status],
          alignSelf: 'center',
          left: '-20px',
          position: 'absolute',
        }} 
      />;
    };

    const nodeLabel = (state: AugmentedStat) => {
      const data = state.data;
    
      return data.friendlyName ?? data.configId;  
    };

    const hasAddButton = (data: PipelineState): data is (PipelineStateSequential<any> | PipelineStateParallel<any>) => 
      (isParallelPipelineState(data) || isSequentialPipelineState(data)) && data.stepTypes.length > 0;

    const treeNodeRef = Vue.ref(null as null | HTMLElement);
    Vue.watch(useElementHover(treeNodeRef), (isHovered) => {
      props.stat.data.isHovered = isHovered;
    });
    
    return () => (
      <div 
        style={{
          display: 'flex', 
          width: '100%',
          height: '30px',
          alignItems: 'center',
          borderBottom: '1px solid var(--steel-2)',
        }}
        ref={treeNodeRef}
        onClick={() => emit('click')}
      >
        { props.callState && progressIcon(callStateToStatus(props.callState)) }
        { props.stat.children.length ? openIcon() : null }
        <span class="mtl-ml text-nowrap text-ellipsis overflow-hidden">{ nodeLabel(props.stat) }</span>
        { props.stat.data.isHovered ? 
          <div class='flex items-center px-2 w-fit justify-end ml-auto'>
            { hasAddButton(props.stat.data) ? 
              <ComboPopup 
                caption={ui.iconFA('plus')}
                items={props.stat.data.stepTypes
                  .map((stepType) => stepType.friendlyName ?? stepType.nqName ?? stepType.configId)
                }
                onSelected={({itemIdx}) => {
                  const data = props.stat.data as PipelineStateSequential<any> | PipelineStateParallel<any>;
                  emit('addNode', {
                    itemId: data.stepTypes[itemIdx].configId,
                    position: data.steps.length,
                  });
                }}
                class='d4-ribbon-item'
              />: null }
            { props.isDraggable ? <IconFA 
              name='grip-vertical' 
              cursor='grab'
              class='d4-ribbon-item'
            />: null }
            { props.isDeletable ? <IconFA 
              name='times' 
              onClick={(e: Event) => {emit('removeNode'); e.stopPropagation();}}
              class='d4-ribbon-item'
            />: null }
          </div> : null }
      </div>
    );    
  },
});

