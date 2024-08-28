import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { defineComponent, h, nextTick, onMounted, onUnmounted, onUpdated, reactive, ref, SlotsType, Teleport, watch } from 'vue';

export const RibbonPanels = defineComponent({
  name: 'RibbonPanels',
  slots: Object as SlotsType<{
    default?: any,
  }>,
  setup(_, {slots}) {
    const elements = reactive([] as HTMLElement[])

    onMounted(async () => {
      await nextTick();

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

    onUnmounted(() => {
      const currentView = grok.shell.v;

      const filteredPanels = currentView
        .getRibbonPanels()
        .filter((panel) => !panel.some((ribbonItem, idx) => ribbonItem === elements[idx].parentElement))

      currentView.setRibbonPanels(filteredPanels);
    });

    return () => 
      slots.default?.().map((slot: any) => <div ref={(el) => addElement(el)}> { slot } </div>) 
    ;
  }
});