import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, KeepAlive, PropType, SlotsType} from 'vue';
import type {DockSpawnTsWebcomponent} from '@datagrok-libraries/webcomponents/src';

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
  setup(_, {slots, emit}) {
    return () => {
      return <dock-spawn-ts 
        style={{'width': '100%'}}
        onPanelClosed={(ev: {detail: any}) => emit('panelClosed', ev.detail)}
      >
        { slots.default?.() }
      </dock-spawn-ts>
    };
  },
});
