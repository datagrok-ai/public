import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { defineComponent, SlotsType } from "vue";
import { TabHeaderStripe } from "./TabHeaderStripe/TabHeaderStripe";
import { TabArea } from "./TabArea/TabArea";
import { HeaderItem } from './shared';

export const Tabs = defineComponent({
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
  slots: Object as SlotsType<{
    default?: any,
    header: HeaderItem
  }>,
  emits: {
    'update:selected': (item: number) => item,
  },
  setup(props, {emit, slots}) {
    console.log(slots.default?.(), ' in tabs');
    return () => (
      <div style={{display: 'flex', flexDirection: 'column'}}>
        <TabHeaderStripe items={props.items} selected={props.selected} onUpdate:selected={(v) => emit('update:selected', v)}>
          {{header: slots.header }}
        </TabHeaderStripe>
        <TabArea selectedIdx={props.selected}> 
          { ...slots.default?.() }
        </TabArea>
      </div>
    )
  }
})