import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {IState} from '@datagrok-libraries/dock-spawn-dg/src/interfaces/IState';
// cannot use DockSpawnTsWebcomponent typings in strict mode
type DockSpawnTsWebcomponent = any

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
    activePanelTitle: Object as Vue.PropType<String>,
  },
  slots: Object as Vue.SlotsType<{
    default?: Vue.VNode[],
  }>,
  emits: {
    'panelClosed': (_element: HTMLElement) => true,
    'panelOpened': (_element: HTMLElement) => true,
    'update:activePanelTitle': (_newPanel: string | null, _prevPanel: string | null) => true,
  },
  methods: {
    'saveLayout': (): string => '',
    'loadLayout': async (_a: string) => void(0),
  },
  setup(props, {slots, emit, expose}) {
    Vue.onRenderTriggered((event) => {
      console.log('DockManager onRenderTriggered', event);
    });

    const dockSpawnRef = Vue.shallowRef<DockSpawnTsWebcomponent | undefined>(undefined);
    const activePanelTitle = Vue.computed(() => props.activePanelTitle);

    const saveLayout = () => dockSpawnRef.value?.saveLayout();
    const loadLayout = (layout: string) => dockSpawnRef.value?.loadLayout(layout);

    expose({saveLayout, loadLayout});

    const onManagerInitFinished = (ev: any) => {
      const manager = ev.detail;
      manager.getElementCallback = async (state: IState) => {
        let aimSlot = null as null | HTMLElement;

        const slots = dockSpawnRef.value.shadowRoot!
          .querySelectorAll(`slot`);
        slots.forEach((slot: any) => {
          const content = (slot.assignedElements() as HTMLElement[])
            .find((el) => el.getAttribute('dock-spawn-title') && el.getAttribute('dock-spawn-title') === state.element);
          if (content) aimSlot = slot;
        });

        return {
          element: aimSlot!,
          title: state.element!,
        };
      };
    }

    return () => (
      <dock-spawn-ts
        style={{'width': '100%'}}
        activePanelTitle={activePanelTitle.value}
        onPanelClosed={(ev: {detail: any}) => emit('panelClosed', ev.detail)}
        onPanelOpened={(ev: {detail: any}) => emit('panelOpened', ev.detail)}
        onActivePanelChanged={(ev: {detail: {newPanel: string | null, prevPanel: string | null}}) => emit('update:activePanelTitle', ev.detail.newPanel, ev.detail.prevPanel)}
        onManagerInitFinished={onManagerInitFinished}
        ref={dockSpawnRef}
      >
        { slots.default?.() }
      </dock-spawn-ts>
    );
  }
});
