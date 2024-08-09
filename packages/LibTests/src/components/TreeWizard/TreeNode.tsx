import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, PropType} from 'vue';
import {AugmentedStat, Status, TreeNodeType} from './types';
import {IconFA} from '@datagrok-libraries/webcomponents-vue/src';
import {OpenIcon} from '@he-tree/vue';

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

const runInParallel = async (children: AugmentedStat[]) => Promise.all(children.map(async (child) => runByTree(child as AugmentedStat)));

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
    parent.status = currentStat.status;
    const parentNextSibling = getNextSibling(parent);
    if (parentNextSibling)
      parentNextSibling.status = `didn't run`;
  }
};

const runByTree = async (currentStat: AugmentedStat): Promise<Status> => {
  const nodeCall = currentStat.data.funcCall;
  currentStat.status = 'running';
  if (nodeCall) {
    try {
      await getCall(nodeCall).call();
      currentStat.status = 'succeeded';
      
      updateNextStepAndParent(currentStat);

      return Promise.resolve(currentStat.status);
    } catch (e) {
      currentStat.status = 'failed';

      return Promise.reject(e);
    }
  } else {  
    return (currentStat.data.order === `parallel` ? 
      runInParallel(currentStat.children): 
      runSequentallly(currentStat.children)
    ).then(() => {
      currentStat.status = 'succeeded';

      updateNextStepAndParent(currentStat);

      return Promise.resolve(currentStat.status);
    }).catch((e) => {
      currentStat.status = 'failed';

      return Promise.reject(e);
    });                      
  }
};

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


export const TreeNode = defineComponent({
  props: {
    stat: {
      type: Object as PropType<AugmentedStat>,
      required: true,
    },
    node: {
      type: Object as PropType<TreeNodeType>,
      required: true,
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
    addNode: (text: String) => text,
    removeNode: () => {},
    click: () => {},
  },
  setup(props, {emit}) {
    if (!props.stat.status)  
      props.stat.status = getInitialStatus(props.stat);
    
    const runIcon = <IconFA 
      name='play'
      onClick={() => runByTree(props.stat)}
      style={{paddingLeft: '4px'}}
    />;
      
    const openIcon = <OpenIcon
      open={props.stat.open}
      class="mtl-mr"
      //@ts-ignore
      onClick={(e) => {props.stat.open = !props.stat.open; e.stopPropagation();}}
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
    
    return () => (
      <div 
        style={{
          display: 'flex', 
          padding: '4px 0px', 
          width: '100%',
        }}
        onMouseover={() => props.stat.isHovered = true} 
        onMouseleave={() => props.stat.isHovered = false} 
        onDragstart={() => props.stat.isHovered = false}
        onClick={() => emit('click')}
        class={props.stat.status === 'locked' ? 'd4-disabled': null}
      >
        { progressIcon(props.stat.status) }
        { props.stat.children.length ? openIcon : null }
        <span class="mtl-ml">{props.node.text}</span>
        { props.isDraggable && props.stat.isHovered && props.stat.status !== 'running' ? <IconFA 
          name='grip-vertical' 
          cursor='grab'
          style={{paddingLeft: '4px'}}
        />: null }
        { props.isDeletable && props.stat.isHovered && props.stat.status !== 'running' ? <IconFA 
          name='times' 
          style={{paddingLeft: '4px'}}
          onClick={(e: Event) => {emit('removeNode'); e.stopPropagation();}}
        />: null }
        { props.isDroppable && props.stat.isHovered && props.stat.status !== 'running' ? <IconFA 
          name='plus' 
          style={{paddingLeft: '4px'}}
          onClick={(e) => {
            emit('addNode', 
              `${props.stat.data.text.includes('phase') ? 'Phase': 'Day'} ${((props.stat.children.length ?? 0) + 1)}`,
            );
            e.stopPropagation();
          }}
        />: null }
        { props.stat.isHovered && props.stat.status !== 'running' ? runIcon: null }
      </div>
    );    
  },
});

