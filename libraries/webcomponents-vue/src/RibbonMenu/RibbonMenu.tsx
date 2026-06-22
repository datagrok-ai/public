import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

const shallowEqualNodes = (a: HTMLElement[], b: HTMLElement[]): boolean =>
  a.length === b.length && a.every((el, i) => el === b[i]);

export const RibbonMenu = Vue.defineComponent({
  name: 'RibbonMenu',
  props: {
    groupName: {
      type: String,
      required: true,
    },
    view: {
      type: DG.ViewBase,
      required: true,
    },
  },
  slots: Object as Vue.SlotsType<{
    default?: Vue.VNode[],
  }>,
  setup(props, {slots}) {
    Vue.onRenderTriggered((event) => {
      console.log('RibbonMenu onRenderTriggered', event);
    });

    const elements = Vue.reactive(new Map<number, HTMLElement>);

    const currentView = Vue.computed(() => Vue.markRaw(props.view));

    let lastApplied: HTMLElement[] = [];
    let lastGroupName: string | undefined = undefined;

    Vue.watch(elements, () => {
      const elementsArray = [...elements.values()];
      // Skip the group rebuild when neither the hosted items nor the group name changed.
      // onBeforeUpdate clears + re-adds the same DOM nodes every re-render; re-pushing is a no-op.
      if (props.groupName === lastGroupName && shallowEqualNodes(elementsArray, lastApplied))
        return;

      currentView.value.ribbonMenu
        .group(props.groupName)
        .clear();

      currentView.value.ribbonMenu
        .group(props.groupName).items(elementsArray, () => {});

      lastApplied = elementsArray;
      lastGroupName = props.groupName;
    }, {flush: 'post'});

    const addElement = (el: Element | null | any, idx: number) => {
      const content = el ? Vue.markRaw(el) : undefined
      if (content)
        elements.set(idx, content);
    };

    Vue.onUnmounted(() => {
      currentView.value.ribbonMenu
        .group(props.groupName)
        .clear();

      currentView.value.ribbonMenu.remove(props.groupName);
    });

    Vue.onBeforeUpdate(() => {
      elements.clear();
    });

    return () =>
      slots.default?.()
        .filter((slot) => slot.type !== Symbol.for('v-cmt'))
        .flatMap((slot) => slot.type === Symbol.for('v-fgt') ? slot.children: slot)
        .map((slot, idx) =>
          <div slot-idx={`${idx}`} style={{width: '100%'}} ref={(el) => addElement(el, idx)}> { slot } </div>,
        );
  },
});
