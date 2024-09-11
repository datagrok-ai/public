import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { defineComponent } from 'vue';
import type {DGMarkdownT} from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-markdown': DGMarkdownT
    }
  }
}

export const MarkDown = defineComponent({
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