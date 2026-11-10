import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useDgView} from '@datagrok-libraries/webcomponents-vue';

type OpenedHistory = {nqName: string, readonly: boolean};

const handedOver = new WeakMap<DG.ViewBase, OpenedHistory>();

// a workflow rebuilds its step view on opening a run; an open Runs panel passes to the next view of the same history
export function useHistoryPanel(nqName: () => string | undefined) {
  const view = useDgView();
  const opened = handedOver.get(view);
  handedOver.delete(view);
  const kept = opened?.nqName === nqName() ? opened : undefined;

  const historyHidden = Vue.ref(!kept);
  const historyReadonly = Vue.ref(!!kept?.readonly);

  Vue.onBeforeUnmount(() => {
    const name = nqName();
    if (!historyHidden.value && name)
      handedOver.set(view, {nqName: name, readonly: historyReadonly.value});
  });

  return {historyHidden, historyReadonly};
}
