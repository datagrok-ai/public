/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Viewer} from '@datagrok-libraries/webcomponents-vue';
import * as Vue from 'vue';

export const ViewerTestApp = Vue.defineComponent({
  name: 'ViewerTestApp',
  setup() {
    const df = Vue.shallowRef<DG.DataFrame | undefined>(undefined);
    const type = Vue.ref<string | undefined>(undefined);
    const viewer = Vue.shallowRef<DG.Viewer | undefined>(undefined);

    let i = 0;
    const datasets = [grok.data.demo.demog(), grok.data.demo.doseResponse(), grok.data.demo.geo()];
    const changeData = () => {
      df.value = Vue.markRaw(datasets[i]);
      i++;
      i%=datasets.length;
    };
    const changeViewer = (v?: DG.Viewer) => {
      viewer.value = v ? Vue.markRaw(v) : undefined;
    };
    let j = 0;
    const types = ['Grid', 'Histogram', 'Line chart', 'Scatter plot'];
    const changeType = () => {
      type.value = types[j];
      j++;
      j%=types.length;
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <button onClick={changeData}>change data</button>
        <button onClick={changeType}>change type</button>
        <Viewer type={type.value} dataFrame={df.value} style={{width: '100%', height: '100%'}} onViewerChanged={changeViewer}></Viewer>
      </div>
    );
  },
});
