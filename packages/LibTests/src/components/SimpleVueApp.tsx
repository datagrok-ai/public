/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
// eslint-disable-next-line
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { Viewer } from '@datagrok-libraries/webcomponents-vue';
import { defineComponent, shallowRef, ref } from 'vue';

export const SimpleTestApp = defineComponent({
  name: 'SimpleTestApp',
  mounted() {
    console.log('SimpleTestApp mounted');
  },
  unmounted() {
    console.log('SimpleTestApp unmounted');
  },
  setup() {
    const df = shallowRef<DG.DataFrame | undefined>(undefined);
    const name = ref<string | undefined>(undefined);
    let i = 0;
    const datasets = [grok.data.demo.demog(), grok.data.demo.doseResponse(), grok.data.demo.geo(), grok.data.demo.geo()];
    function changeData() {
      df.value = datasets[i];
      i++;
      i%=datasets.length;
    }
    let j = 0;
    const types = ['Grid', 'Histogram', 'Line chart', 'Scatter plot'];
    function changeType() {
      name.value = types[j];
      j++;
      j%=types.length;
    }
    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <button onClick={() => changeData()}>change data</button>
        <button onClick={() => changeType()}>change type</button>
        <Viewer name={name.value} value={df.value}></Viewer>
      </div>
    );
  }
});
