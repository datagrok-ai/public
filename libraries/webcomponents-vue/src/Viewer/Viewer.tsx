import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {ViewerT} from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const Viewer = Vue.defineComponent({
  name: 'Viewer',
  props: {
    type: String,
    dataFrame: DG.DataFrame,
    options: Object as Vue.PropType<Record<string, string | boolean | string[]>>,
  },
  emits: {
    viewerChanged: (_v: DG.Viewer<any> | undefined) => true,
  },
  setup(props, {emit}) {
    Vue.onRenderTriggered((event) => {
      console.log('Viewer onRenderTriggered', event);
    });

    const currentDf = Vue.computed(() => props.dataFrame ? Vue.markRaw(props.dataFrame) : undefined);
    const options = Vue.computed(() => props.options ? Vue.markRaw(props.options) : undefined);
    const type = Vue.computed(() => props.type);
    const viewerChangedCb = (event: any) => {
      emit('viewerChanged', event.detail ? Vue.markRaw(event.detail) : undefined);
    };
    return () => <Vue.KeepAlive>
      <dg-viewer
        type={type.value}
        options={options.value}
        dataFrame={currentDf.value}
        onViewerChanged={viewerChangedCb}
        style={{display: 'block', flexGrow: '1'}}
      >
      </dg-viewer>
    </Vue.KeepAlive>;
  },
});
