import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {IconFA, SplitH,  Tabs} from '@datagrok-libraries/webcomponents-vue';
import * as Vue from 'vue';

export const ElementsApp = Vue.defineComponent({
  name: 'ElementsApp',
  setup() {
    const resize = Vue.ref(true);

    const items = Vue.ref([{label: 'Tab 1'}, {label: 'Tab 2'}, {label: 'Tab 3'}]);
    const selected = Vue.ref(1);

    return () => (
      <keep-alive>
        <div style={{width: '100%', height: '100%'}}>
          <Tabs items={items.value} selected={selected.value} onUpdate:selected={(v) => selected.value = v}>
            {{
              header: (item: any) => <span><IconFA name='info'/> {item.label} </span>,
              default: () => [
                (<div> <span> 1 </span> </div>),
                (<div> <span> 2 </span> </div>),
              ],
            }}
          </Tabs>

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
