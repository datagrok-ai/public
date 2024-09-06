import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, KeepAlive, onMounted, PropType, ref, SlotsType} from 'vue';
import {DockSpawnTsWebcomponent} from '@datagrok-libraries/webcomponents/src';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dock-spawn-ts': DockSpawnTsWebcomponent
    }
  }
}

export const DockManager = defineComponent({
  name: 'DockManager',
  slots: Object as SlotsType<{
    default?: any,
  }>,
  emits: {
    panelClosed: (element: HTMLElement) => element,
  },
  methods: {
    saveLayout: () => {},
    loadLayout: () => {}
  },
  setup(_, {slots, emit, expose}) {
    let dockSpawnRef = ref(null as DockSpawnTsWebcomponent | null);
    const saveLayout = () => {
      if (!dockSpawnRef.value) return;

      localStorage.setItem('test-layout', dockSpawnRef.value.saveLayout())
    }

    const loadLayout = () => {
      const savedLayout = localStorage.getItem('test-layout');

      if (!dockSpawnRef.value || !savedLayout) return;

      dockSpawnRef.value.loadLayout(savedLayout);
    }
    expose({
      'saveLayout': saveLayout,
      'loadLayout': loadLayout
    })

    return () => {
      return <dock-spawn-ts 
        style={{'width': '100%'}}
        onPanelClosed={(ev: {detail: any}) => emit('panelClosed', ev.detail)}
        ref={dockSpawnRef}
      >
        { slots.default?.() }
      </dock-spawn-ts>
    };
  },
});
