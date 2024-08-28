import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { defineComponent, nextTick, onMounted, onUpdated, reactive, ref, SlotsType, Teleport, watch } from 'vue';

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

    return () => 
      slots.default?.().map((slot: any) => <div ref={(el) => addElement(el)}> { slot } </div>) 
    ;
  }
});