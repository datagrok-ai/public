import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export const TabArea = Vue.defineComponent({
  props: {
    selectedIdx: {
      type: Number,
      default: 0
    }
  },
  setup(props, {slots}){
    return () => (
      <div class={'d4-tab-content ui-box'}>
        { slots.default?.().at(props.selectedIdx) ?? null }
      </div>
    )
  }
});