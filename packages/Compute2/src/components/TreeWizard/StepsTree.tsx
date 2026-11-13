import * as Vue from 'vue';
import {Draggable, dragContext} from '@he-tree/vue';
import '@he-tree/vue/style/default.css';
import '@he-tree/vue/style/material-design.css';
import {
  isDynamicPipelineState, PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {AugmentedStat} from './types';
import {dropIndex, focusCandidates, isDeletable, isEachDraggable, stepMoveIndex} from './tree-utils';
import {isKey, TAB_CYCLE_CLASS} from '../../composables/use-keyboard';

/**
 * The workflow's steps on he-tree. Drops and Alt+Up/Down move a step under the drag rules and Delete removes one, all
 * through the parent; the collapsed nodes and the focused row carry over he-tree's rebuilds. The rows come from the
 * default slot.
 */
export const StepsTree = Vue.defineComponent({
  name: 'StepsTree',
  props: {
    treeState: {
      type: Object as Vue.PropType<PipelineState>,
      required: true,
    },
    // the collapsed nodes are kept by the parent, as the tree unmounts while its panel is hidden
    closed: {
      type: Object as Vue.PropType<Set<string>>,
      required: true,
    },
  },
  emits: {
    choose: (_state: PipelineState) => true,
    move: (_uuid: string, _index: number) => true,
    remove: (_state: PipelineState) => true,
  },
  slots: Object as Vue.SlotsType<{default: {stat: AugmentedStat, isDraggable: boolean}}>,
  setup(props, {emit, slots, expose}) {
    const tree = Vue.shallowRef(null as InstanceType<typeof Draggable> | null);

    const statOf = (el: Element | null) => {
      const row = el?.closest<HTMLElement>('.tree-node');
      return row && tree.value?.$el.contains(row) ? tree.value.getNodeByElement(row) as AugmentedStat | null : null;
    };
    expose({stepAt: (el: Element | null) => statOf(el)?.data});

    // he-tree rebuilds its rows from every new state, so the closed nodes and the focused row carry over by uuid
    let oldFocus: string[] = [];

    Vue.watch(() => props.treeState, () => {
      if (!tree.value)
        return;
      props.closed.clear();
      for (const stat of tree.value.statsFlat as AugmentedStat[]) {
        if (!stat.open)
          props.closed.add(stat.data.uuid);
      }
      const focused = statOf(document.activeElement);
      oldFocus = focused ? focusCandidates(focused) : [];
    });

    // a row moved, removed or added during the rebuild loses the focus, which goes back to the closest remaining step
    Vue.watch(() => props.treeState, () => {
      const candidates = oldFocus;
      oldFocus = [];
      const el = tree.value?.$el as HTMLElement | undefined;
      if (!candidates.length || !el || el.contains(document.activeElement))
        return;
      const rows = new Map(([...el.querySelectorAll('.tree-node')] as HTMLElement[])
        .map((row) => [tree.value!.getNodeByElement(row)?.data.uuid, row]));
      candidates.map((uuid) => rows.get(uuid)).find((row) => !!row)?.focus();
    }, {flush: 'post'});

    const restoreClosed = (stat: AugmentedStat) => {
      if (props.closed.has(stat.data.uuid))
        stat.open = false;
      return stat;
    };

    const isEachDroppable = (stat: AugmentedStat) => {
      const draggedStep = dragContext?.startInfo?.dragNode as AugmentedStat | undefined;
      return isDynamicPipelineState(stat.data) &&
        (!draggedStep || stat.data.uuid === draggedStep.parent?.data.uuid);
    };

    const onAfterDrop = () => {
      const draggedStep = dragContext.startInfo?.dragNode as AugmentedStat | undefined;
      if (draggedStep) {
        const oldIndex = dragContext.startInfo.indexBeforeDrop;
        const newIndex = dropIndex(oldIndex, dragContext.targetInfo.indexBeforeDrop);
        if (oldIndex !== newIndex)
          emit('move', draggedStep.data.uuid, newIndex);
      }
      // he-tree skips inbound modelValue rebuilds while its `dragNode` data is set. A slow
      // frame or a debugger pause during drag-end can strand `dragNode` truthy, which then
      // freezes every later structural update (added/removed steps stop rendering). The
      // library clears it in its own `dragend` handler, but that path is unreliable under
      // interruption; `after-drop` fires reliably, so clear it here too.
      const inst = tree.value as any;
      if (inst) {
        inst.dragNode = null;
        inst.dragOvering = false;
      }
    };

    // he-tree's own Alt+arrow moves ignore the drag rules and never reach the driver, so the tree takes these keys
    // first: Alt+Up/Down moves a step the way a drop does, Delete removes it
    const onKeydown = (e: KeyboardEvent) => {
      const stat = statOf(e.target as Element);
      if (!stat)
        return;
      if (e.altKey && e.key.startsWith('Arrow')) {
        e.preventDefault();
        e.stopPropagation();
        const index = e.key === 'ArrowUp' || e.key === 'ArrowDown' ?
          stepMoveIndex(stat, e.key === 'ArrowUp' ? -1 : 1) : undefined;
        if (index !== undefined)
          emit('move', stat.data.uuid, index);
      } else if (isKey(e, ['Delete']) && isDeletable(stat)) {
        e.preventDefault();
        emit('remove', stat.data);
      }
    };

    return () => (
      <Draggable
        class={`ui-div mtl-tree p-2 overflow-scroll h-full ${TAB_CYCLE_CLASS}`}
        style={{paddingLeft: '25px'}}
        rootDroppable={false}
        treeLine
        updateBehavior='disabled'
        childrenKey='steps'
        nodeKey={(stat: AugmentedStat) => stat.data.uuid}
        statHandler={restoreClosed}
        ref={tree}
        modelValue={[props.treeState]}
        eachDraggable={isEachDraggable}
        eachDroppable={isEachDroppable}
        onAfter-drop={onAfterDrop}
        {...{onKeydownCapture: onKeydown}}
        onClick:node={(stat: AugmentedStat) => emit('choose', stat.data)}
        onClose:node={(stat: AugmentedStat) => props.closed.add(stat.data.uuid)}
        onOpen:node={(stat: AugmentedStat) => props.closed.delete(stat.data.uuid)}
      >
        {({stat}: {stat: AugmentedStat}) => slots.default?.({stat, isDraggable: !!tree.value?.isDraggable(stat)})}
      </Draggable>
    );
  },
});
