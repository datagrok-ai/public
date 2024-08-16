import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { defineComponent, PropType, SlotsType } from 'vue';

type HeaderItem  = {label: string, key?: string};

export const TabHeaderStripe = defineComponent({
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
    header: HeaderItem
  }>,
  emits: {
    'update:selected': (item: number) => item,
  },
  setup(props, {slots, emit}) {
    return () => (
      <div class={'d4-tab-header-stripe'}> 
      { 
        props.items.map((item, idx) => (
          <div 
            class={{'d4-tab-header': true, 'selected': (idx === props.selected)}} 
            key={item.key ?? item.label} 
            onClick={() => emit('update:selected', idx)}
          >
            { slots.header?.(item) ?? item.label }
          </div>
        ))
      }
      </div>
    )
  }
})