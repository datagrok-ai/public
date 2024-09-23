import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {DGMarkdownT} from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-markdown': DGMarkdownT
    }
  }
}

export const MarkDown = Vue.defineComponent({
  name: 'MarkDown',
  props: {
    markdown: {
      type: String,
      required: true,
    }
  },
  setup(props) {
    return () => <dg-markdown content={props.markdown}/>
  }
});