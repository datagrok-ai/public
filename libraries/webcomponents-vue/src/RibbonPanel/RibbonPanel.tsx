import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const RibbonPanel = Vue.defineComponent({
  name: 'RibbonPanel',
  props: {
    view: {
      type: DG.ViewBase,
      required: true,
    },
  },
  slots: Object as Vue.SlotsType<{
    default?: any,
  }>,
  setup(props, {slots}) {
    const elements = Vue.reactive(new Map<number, HTMLElement>);

    const currentView = Vue.shallowRef(props.view);

    Vue.watch(elements, () => {
      const elementsArray = [...elements.values()];
      const filteredPanels = currentView.value
        .getRibbonPanels()
        .filter((panel) => !panel.some((ribbonItem) => elementsArray.includes(ribbonItem.children[0] as HTMLElement)));
      currentView.value.setRibbonPanels(filteredPanels);

      currentView.value.setRibbonPanels([
        currentView.value.getRibbonPanels().flat(),
        elementsArray,
      ]);

      // Workaround for ui.comboPopup elements. It doesn't work if it is not a direct child of '.d4-ribbon-item'
      elementsArray.forEach((elem) => {
        const content = ((elem.firstChild?.nodeType !== Node.TEXT_NODE) ? elem.firstChild: elem.firstChild.nextSibling) as HTMLElement | null;

        if (content && content.tagName.toLowerCase().includes('dg-combo-popup')) {
          content.classList.add('d4-ribbon-item');
          elem.parentElement?.classList.remove('d4-ribbon-item');
        }
      });
    });

    const addElement = (el: Element | null | any, idx: number) => {
      const content = el;
      if (content)
        elements.set(idx, content);
    };

    Vue.onUnmounted(async () => {
      const elementsArray = [...elements.values()];
      const filteredPanels = currentView.value
        .getRibbonPanels()
        .filter((panel) => !panel.some((ribbonItem) => elementsArray.includes(ribbonItem.children[0] as HTMLElement)));

      await Vue.nextTick();
      if (!currentView.value.closing)
        currentView.value.setRibbonPanels(filteredPanels);
    });

    return () =>
      slots.default?.().filter((slot: any) => slot.type !== Symbol.for('v-cmt')).map((slot: any, idx: number) =>
        <div slot-idx={`${idx}`} ref={(el) => addElement(el, idx)}> { slot } </div>,
      );
  },
});
