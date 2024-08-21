
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {defineComponent} from 'vue';
import {History} from './History';

export const HistoryApp = defineComponent({
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
