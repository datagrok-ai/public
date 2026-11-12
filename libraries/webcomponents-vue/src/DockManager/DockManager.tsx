import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

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
    activePanelTitle: String,
    preferredPanelTitle: String,
  },
  slots: Object as Vue.SlotsType<{
    default?: Vue.VNode[],
  }>,
  emits: {
    'panelClosed': (_element: HTMLElement) => true,
    'panelOpened': (_element: HTMLElement) => true,
    'update:activePanelTitle': (_newPanel: string | null, _prevPanel: string | null) => true,
    'tabClicked': (_title: string | null) => true,
  },
  setup(props, {slots, emit}) {
    const activePanelTitle = Vue.computed(() => props.activePanelTitle);
    const preferredPanelTitle = Vue.computed(() => props.preferredPanelTitle);

    return () => (
      <dock-spawn-ts
        style={{'width': '100%'}}
        activePanelTitle={activePanelTitle.value}
        preferredPanelTitle={preferredPanelTitle.value}
        onPanelClosed={(ev: {detail: any}) => emit('panelClosed', ev.detail)}
        onPanelOpened={(ev: {detail: any}) => emit('panelOpened', ev.detail)}
        onActivePanelChanged={(ev: {detail: {newPanel: string | null, prevPanel: string | null}}) =>
          emit('update:activePanelTitle', ev.detail.newPanel, ev.detail.prevPanel)}
        onTabClicked={(ev: {detail: {title: string | null}}) => emit('tabClicked', ev.detail.title)}
      >
        { slots.default?.() }
      </dock-spawn-ts>
    );
  },
});
