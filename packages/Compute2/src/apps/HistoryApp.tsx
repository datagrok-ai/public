
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {History} from '../components/History/History';

export const HistoryApp = Vue.defineComponent({
  name: 'HistoryApp',
  setup: () => {
    
    return () => <History
      func={DG.Func.byName('Compute:ObjectCooling')}
      showActions={true}
      showBatchActions={true}
      isHistory={true}
    />;
  },
});
