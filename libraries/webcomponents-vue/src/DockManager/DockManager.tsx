import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {DockSpawnTsWebcomponent} from '@datagrok-libraries/webcomponents';
import { IState } from '@datagrok-libraries/webcomponents/vendor/dock-spawn-ts/lib/interfaces/IState';
import { whenever } from '@vueuse/core';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dock-spawn-ts': DockSpawnTsWebcomponent
    }
  }
}

export const DockManager = Vue.defineComponent({
  name: 'DockManager',
  props: {
    layoutStorageName: Object as Vue.PropType<String>
  },
  slots: Object as Vue.SlotsType<{
    default?: Vue.VNode[],
  }>,
  emits: {
    panelClosed: (element: HTMLElement) => element,
    activePanelChanged: (panelTitle: string) => panelTitle
  },
  methods: {
    saveLayout: () => {},
    loadLayout: () => {}
  },
  setup(props, {slots, emit, expose}) {
    let dockSpawnRef = Vue.ref(null as DockSpawnTsWebcomponent | null);

    whenever(dockSpawnRef, async () => {
      if (!dockSpawnRef.value?.dockManager) await Vue.nextTick();

      dockSpawnRef.value!.dockManager.getElementCallback = async (state: IState) => {
        let aimSlot = null as null | HTMLElement;

        const slots = dockSpawnRef.value!.shadowRoot!
          .querySelectorAll(`slot`);
        slots.forEach((slot) => {
          const content = (slot.assignedElements() as HTMLElement[])
            .find((el) => el.getAttribute('dock-spawn-title') && el.getAttribute('dock-spawn-title') === state.element); 
          if (content) aimSlot = slot;
        })

        return { 
          element: aimSlot!,
          title: state.element!,
        }
      };
    }, {once: true})

    const saveLayout = () => {
      if (!dockSpawnRef.value || !props.layoutStorageName) return;

      localStorage.setItem(props.layoutStorageName.toString(), dockSpawnRef.value.saveLayout())
    }

    const loadLayout = async () => {
      if (!dockSpawnRef.value || !props.layoutStorageName) return;

      const savedLayout = localStorage.getItem(props.layoutStorageName.toString());

      if (!savedLayout) return;

      await dockSpawnRef.value.loadLayout(savedLayout)
    }
    expose({
      'saveLayout': saveLayout,
      'loadLayout': loadLayout
    })

    return () => {
      return <dock-spawn-ts 
        style={{'width': '100%'}}
        onPanelClosed={(ev: {detail: any}) => emit('panelClosed', ev.detail)}
        onActivePanelChanged={(ev: {detail: any}) => emit('activePanelChanged', ev.detail)}
        ref={dockSpawnRef}
      >
        { slots.default?.() }
      </dock-spawn-ts>
    };
  },
});
