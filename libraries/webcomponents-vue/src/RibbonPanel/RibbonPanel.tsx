import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const RibbonPanel = Vue.defineComponent({
  name: 'RibbonPanel',
  slots: Object as Vue.SlotsType<{
    default?: any,
  }>,
  setup(_, {slots}) {
    const elements = Vue.reactive([] as HTMLElement[])

    Vue.onMounted(async () => {
      await Vue.nextTick();

      const currentView = grok.shell.v;
      currentView.setRibbonPanels([
        currentView.getRibbonPanels().flat(),
        elements,
      ]);
    })    

    const addElement = (el: Element | null | any) => {
      if (el)
        elements.push(el)
    }

    Vue.onUnmounted(() => {
      const currentView = grok.shell.v;

      const filteredPanels = currentView
        .getRibbonPanels()
        .filter((panel) => !panel.some((ribbonItem, idx) => ribbonItem === elements[idx].parentElement))

      currentView.setRibbonPanels(filteredPanels);
    });

    return () => 
      slots.default?.().map((slot: any) => <div ref={(el) => addElement(el)}> { slot } </div>) 
  }
});