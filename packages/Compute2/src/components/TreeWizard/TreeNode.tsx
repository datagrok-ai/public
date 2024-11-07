import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {AugmentedStat, Status} from './types';
import {ComboPopup, IconFA} from '@datagrok-libraries/webcomponents-vue';
import {OpenIcon} from '@he-tree/vue';
import {
  isParallelPipelineState,
  isSequentialPipelineState, PipelineState,
  PipelineStateParallel, PipelineStateSequential,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {useElementHover} from '@vueuse/core';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';

const statusToIcon = {
  ['next']: 'arrow-right',
  [`pending`]: 'circle',
  ['running']: 'hourglass-half',
  ['succeeded']: 'check-circle',
  ['partially succeeded']: 'dot-circle',
  ['failed']: 'times-circle',
} as Record<Status, string>;

const statusToColor = {
  ['next']: 'blue',
  [`pending`]: 'gray',
  ['running']: 'blue',
  ['succeeded']: 'green',
  ['partially succeeded']: 'orange',
  ['failed']: 'red',
} as Record<Status, string>;

const statusToTooltip = {
  [`next`]: `This step is avaliable to run`,
  ['pending']: 'This step has pending dependencies',
  ['running']: 'This step is running',
  ['succeeded']: 'This step is succeeded',
  ['partially succeeded']: 'This step has issues',
  ['failed']: 'Run failed',
} as Record<Status, string>;

const statesToStatus = (
  callState: FuncCallStateInfo,
  validationsState?: Record<string, ValidationResult>,
  consistencyStates?: Record<string, ConsistencyInfo>,
): Status => {
  if (callState.isRunning) return 'running';
  if (callState.pendingDependencies?.length) return 'pending';
  if (!callState.isOutputOutdated) {
    return getCallSuccess(validationsState, consistencyStates);
  } else {
    if (callState.runError) return 'failed';
  };
  return 'next';
};

const getToolTip = (status: Status, isReadonly: boolean) => {
  if (!isReadonly || status !== 'next') return statusToTooltip[status];
  return 'This step is locked';
}

const getCallSuccess = (validationsState?: Record<string, ValidationResult>, consistencyStates?: Record<string, ConsistencyInfo>) => {
  const firstError = Object.values(validationsState || {}).find(val => (val.errors?.length || val.warnings?.length))
  const firstInconsistency = Object.values(consistencyStates || {}).find(val => val.inconsistent);
  if (firstError || firstInconsistency)
    return 'partially succeeded';
  return 'succeeded';
}

export const TreeNode = Vue.defineComponent({
  props: {
    stat: {
      type: Object as Vue.PropType<AugmentedStat>,
      required: true,
    },
    callState: {
      type: Object as Vue.PropType<FuncCallStateInfo>,
    },
    validationStates: {
      type: Object as Vue.PropType<Record<string, ValidationResult>>,
    },
    consistencyStates: {
      type: Object as Vue.PropType<Record<string, ConsistencyInfo>>,
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
    isReadonly: {
      type: Boolean,
    }
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

    const progressIcon = (status: Status, isReadOnly: boolean) => {
      return <IconFA
        name={isReadOnly ? 'lock' : statusToIcon[status]}
        animation={status === `running` ? 'spin': null}
        tooltip={getToolTip(status, isReadOnly) ?? null}
        style={{
          color: statusToColor[status],
          alignSelf: 'center',
          left: '-20px',
          position: 'absolute',
        }}
      />;
    };

    const nodeLabel = (state: AugmentedStat) =>
      state.data.friendlyName ?? state.data.configId;

    const hasAddButton = (data: PipelineState): data is (PipelineStateSequential<any> | PipelineStateParallel<any>) =>
      (isParallelPipelineState(data) || isSequentialPipelineState(data)) && data.stepTypes.length > 0 && !data.isReadonly;

    const treeNodeRef = Vue.ref(null as null | HTMLElement);
    const isHovered = useElementHover(treeNodeRef);

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
        { props.callState && progressIcon(statesToStatus(props.callState, props.validationStates, props.consistencyStates), props.isReadonly) }
        { props.stat.children.length ? openIcon() : null }
        <span class="mtl-ml text-nowrap text-ellipsis overflow-hidden">{ nodeLabel(props.stat) }</span>
        { (isHovered.value) ?
          <div class='flex items-center px-2 w-fit justify-end ml-auto'>
            { hasAddButton(props.stat.data) ?
              <ComboPopup
                caption={ui.iconFA('plus')}
                items={props.stat.data.stepTypes
                  .map((stepType) => stepType.friendlyName || stepType.nqName || stepType.configId)
                }
                onSelected={({itemIdx}) => {
                  const data = props.stat.data as PipelineStateSequential<any> | PipelineStateParallel<any>;
                  emit('addNode', {
                    itemId: data.stepTypes[itemIdx].configId,
                    position: data.steps.length,
                  });
                  isHovered.value = false;
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
