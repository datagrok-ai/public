import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BigButton, Button, IconFA, SplitH, TabHeaderStripe} from '@datagrok-libraries/webcomponents-vue/src';
import {defineComponent, ref} from 'vue';

export const VueElementsTestApp = defineComponent({
  name: 'VueElementsTestApp',
  setup() {
    const resize = ref(true);

    const items = ref([{label: 'Tab 1'}, {label: 'Tab 2'}, {label: 'Tab 3'}]);
    const selected = ref(1);

    return () => (
      <keep-alive>
        <div style={{width: '100%', height: '100%'}}>
          <TabHeaderStripe items={items.value} selected={selected.value} onUpdate:selected={(v) => selected.value = v}/>
        
          <SplitH resize={resize.value}>
            <button onClick={() => items.value = [{label: 'Test anme'}, {label: 'Seocnd tab'}]} is='dg-button'>
            Change tab name
            </button>
            <button onClick={() => {
              resize.value = !resize.value;
            }} is='dg-big-button'>click me</button>  
          </SplitH>  
        </div>
      </keep-alive>
    );
  },
});
