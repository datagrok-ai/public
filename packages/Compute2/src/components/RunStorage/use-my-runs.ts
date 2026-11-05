import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {merge, from, of} from 'rxjs';
import {catchError, switchMap} from 'rxjs/operators';
import {from as fromRef, useSubscription} from '@vueuse/rxjs';
import {listMyModelRuns, runsChanged$} from './run-storage';

/** My runs of the models, loaded again when the models change and when runs change; the last load wins. */
export function useMyRuns(models: Vue.MaybeRefOrGetter<DG.Func[]>, stepHistory: Vue.MaybeRefOrGetter<boolean>) {
  const runs = Vue.shallowRef<DG.FuncCall[]>([]);
  const count = Vue.ref<number>();
  const isLoading = Vue.ref(false);
  const error = Vue.ref<string>();

  const key = Vue.computed(() =>
    `${Vue.toValue(models).map((func) => func.nqName).join(',')}|${Vue.toValue(stepHistory)}`);

  // resolves undefined when there is nothing to list or the load failed
  const load = () => {
    const funcs = Vue.toValue(models);
    isLoading.value = funcs.length > 0;
    count.value = undefined;
    error.value = undefined;
    if (funcs.length === 0)
      return of(undefined);
    return from(listMyModelRuns(funcs, Vue.toValue(stepHistory))).pipe(catchError((e) => {
      error.value = String(e?.message ?? e);
      return of(undefined);
    }));
  };

  useSubscription(merge(fromRef(key, {immediate: true}), runsChanged$).pipe(switchMap(load)).subscribe((loaded) => {
    runs.value = (loaded ?? []).map((run) => Vue.markRaw(run));
    count.value = loaded?.length;
    isLoading.value = false;
  }));

  return {
    runs: runs as Readonly<typeof runs>,
    count: count as Readonly<typeof count>,
    isLoading: isLoading as Readonly<typeof isLoading>,
    error: error as Readonly<typeof error>,
  };
}
