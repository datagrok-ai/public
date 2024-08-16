import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { defineComponent } from 'vue';

export const TabArea = defineComponent({
  props: {
    selectedIdx: {
      type: Number,
      default: 0
    }
  },
  setup(props, {slots}){
    console.log(slots.default?.(), ' in tabArea');
    return () => (
      <div class={'d4-tab-content ui-box'}>
        { slots.default?.().at(props.selectedIdx) ?? null }
      </div>
    )
  }
});