import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {v4 as uuidv4} from 'uuid';

const shallowEqualNodes = (a: HTMLElement[], b: HTMLElement[]): boolean =>
  a.length === b.length && a.every((el, i) => el === b[i]);

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
    Vue.onRenderTriggered((event) => {
      console.log('RibbonPanel onRenderTriggered', event);
    });

    const elements = Vue.reactive(new Map<number, HTMLElement>);
    const currentView = Vue.computed(() => Vue.markRaw(props.view));
    const uuid = uuidv4();

    const isInstanceManagedPanel = (panel: HTMLElement[]) => {
      return panel.some((ribbonItem) => (ribbonItem.children[0] as any)?.RibbonPanelUUID === uuid);
    };

    const filterEmptyPanels = (panels: HTMLElement[][]) => {
      return panels.filter((panel) => panel.some((ribbonItem) => (ribbonItem.children?.length)));
    };

    Vue.watch(elements, () => {
      const panels = currentView.value.getRibbonPanels();
      const existingPanelIdx = panels.findIndex(isInstanceManagedPanel);
      const panel = [...elements.values()];

      // Skip the DG round-trip when our panel is already present with the same nodes.
      // Re-renders re-run the ref callbacks with identical DOM nodes; re-pushing them is a no-op.
      if (existingPanelIdx >= 0 && shallowEqualNodes(panels[existingPanelIdx], panel))
        return;

      // Nothing to add/remove
      if (existingPanelIdx < 0 && panel.length == 0) {
        currentView.value.setRibbonPanels(filterEmptyPanels(panels));
        return;
      }

      if (existingPanelIdx >= 0 && panel.length == 0) // Remove panel
        panels.splice(existingPanelIdx, 1);
      else if (existingPanelIdx >= 0) // Replace panel
        panels.splice(existingPanelIdx, 1, panel);
      else // Add new panel
        panels.push(panel);

      currentView.value.setRibbonPanels(filterEmptyPanels(panels));
    });

    const addElement = (el: Element | null | any, idx: number) => {
      const content = el ? Vue.markRaw(el) : undefined;
      if (content) {
        (content).RibbonPanelUUID = uuid;
        elements.set(idx, content);
      } else
        elements.delete(idx);
    };

    Vue.onUnmounted(async () => {
      const notInstanceManagedPanels = currentView.value
        .getRibbonPanels()
        .filter((panel) => !isInstanceManagedPanel(panel));

      await Vue.nextTick();
      if (!currentView.value.closing)
        currentView.value.setRibbonPanels(filterEmptyPanels(notInstanceManagedPanels));
    });

    return () =>
      slots.default?.().filter((slot: any) => slot.type !== Symbol.for('v-cmt')).map((slot: any, idx: number) =>
        <div slot-idx={`${idx}`} ref={(el) => addElement(el, idx)}> { slot } </div>,
      );
  },
});
