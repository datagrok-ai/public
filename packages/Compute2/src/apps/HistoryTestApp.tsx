
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {History} from '../components/History/History';

export const HistoryTestApp = Vue.defineComponent({
  name: 'HistoryTestApp',
  setup: () => {
    return () => <History
      func={DG.Func.byName('Compute2:ObjectCooling2')}
      showActions={true}
      showBatchActions={true}
      isHistory={true}
    />;
  },
});
