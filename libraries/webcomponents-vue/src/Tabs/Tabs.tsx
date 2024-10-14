import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import { TabHeaderStripe } from "./TabHeaderStripe/TabHeaderStripe";
import { TabArea } from "./TabArea/TabArea";
import { HeaderItem } from './shared';

export const Tabs = Vue.defineComponent({
  props: {
    items: {
      type: Array<HeaderItem>,
      required: true,
    },
    selected: {
      type: Number,
      default: 0
    },
  },
  slots: Object as Vue.SlotsType<{
    default: any,
    header: HeaderItem,
    stripePrefix: any,
    stripePostfix: any,
  }>,
  emits: {
    'update:selected': (item: number) => item,
  },
  setup(props, {emit, slots}) {
    return () => (
      <div style={{display: 'flex', flexDirection: 'column', height: '100%'}}>
        <TabHeaderStripe items={props.items} selected={props.selected} onUpdate:selected={(v) => emit('update:selected', v)}>
          {{
            header: slots.header,
            prefix: slots.stripePrefix,
            postfix: slots.stripePostfix,
          }}
        </TabHeaderStripe>
        <TabArea selectedIdx={props.selected}> 
          { ...slots.default?.() }
        </TabArea>
      </div>
    )
  }
})