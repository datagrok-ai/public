import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const RibbonMenu = Vue.defineComponent({
  name: 'RibbonMenu',
  props: {
    groupName: {
      type: String,
      required: true,
    },
    // isChecked: {
    //   type: Object as Vue.PropType<(item: HTMLElement) => boolean>,
    // },
  },
  slots: Object as Vue.SlotsType<{
    default?: any,
  }>,
  setup(props, {slots}) {
    const elements = Vue.reactive(new Map<number, HTMLElement>)

    Vue.watch(elements, async () => {
      await Vue.nextTick();

      const currentView = grok.shell.v;

      const elementsArray = [...elements.values()];
      currentView.ribbonMenu
        .group(props.groupName)
        .clear();

      currentView.ribbonMenu
        .group(props.groupName).items(elementsArray, () => {});
    })    

    const addElement = (el: Element | null | any, idx: number) => {
      const content = el;
      if (content)
        elements.set(idx, content);
    }

    Vue.onUnmounted(() => {
      const currentView = grok.shell.v;

      currentView.ribbonMenu
        .group(props.groupName)
        .clear();

      currentView.ribbonMenu.remove(props.groupName)
    });

    return () => 
      slots.default?.().filter((slot: any) => slot.type !== Symbol.for('v-cmt')).map((slot: any, idx: number) => 
        <div slot-idx={`${idx}`} ref={(el) => addElement(el, idx)}> { slot } </div>
      ) 
  }
});