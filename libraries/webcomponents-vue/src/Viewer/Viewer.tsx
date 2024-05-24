import * as grok from 'datagrok-api/grok';
// eslint-disable-next-line
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { defineComponent, ref } from 'vue';
import type { ViewerT } from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const Viewer = defineComponent({
  name: 'Viewer',
  props: {
    name: String,
    value: DG.DataFrame,
    viewer: DG.Viewer,
  },
  setup(props) {
    return () => {
      const viewer = <dg-viewer name={props.name} value={props.value} viewer={props.viewer}></dg-viewer>;
      return (
        <keep-alive>
          { viewer }
        </keep-alive>
      );
    };
  }
});
