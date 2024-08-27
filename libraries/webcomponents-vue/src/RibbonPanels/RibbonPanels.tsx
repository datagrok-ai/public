import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { defineComponent, nextTick, onMounted, onUpdated, ref, SlotsType, Teleport, watch } from 'vue';

export const RibbonPanels = defineComponent({
  name: 'RibbonPanels',
  slots: Object as SlotsType<{
    default?: any,
  }>,
  setup(_, {slots}) {
    const elements = ref(null as null | HTMLElement)

    watch(elements, (newVal) => {
      if (!newVal) return;

      const currentView = grok.shell.v;
      currentView.setRibbonPanels([
        currentView.getRibbonPanels().flat(),
        Array.from(newVal.children) as HTMLElement[],
      ])
    })

    return () => <div ref={elements}>
      { slots.default?.() }
    </div>;
  }
});