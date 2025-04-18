import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {History} from '../components/History/History';
import {BehaviorSubject, Subject} from 'rxjs';
import {useObservable} from '@vueuse/rxjs';

export const HistoryApp = Vue.defineComponent({
  name: 'HistoryApp',
  props: {
    name: {
      type: String,
      required: true,
    },
    showHistory: {
      type: BehaviorSubject<Boolean>,
      required: true,
    },
    // cannot listen to emited vue events outside of vue
    updateFCBus: {
      type: Subject<DG.FuncCall>,
      required: true,
    },
  },
  setup(props) {
    const showHistory = useObservable(props.showHistory);
    const func = Vue.computed(() => {
      const f = DG.Func.byName(props.name);
      return f ? Vue.markRaw(f) : undefined;
    })
    return () => (
      (showHistory.value && func.value) ?
        <History
          func={func.value}
          showActions={true}
          showBatchActions={true}
          isHistory={true}
          onRunChosen={(fc) => props.updateFCBus.next(Vue.markRaw(fc))}
        /> :
        <div></div>
    );
  },
});
