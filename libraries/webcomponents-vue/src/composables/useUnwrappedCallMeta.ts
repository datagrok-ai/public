import * as Vue from 'vue';
import {BehaviorSubject, merge} from 'rxjs';
import {map, tap} from 'rxjs/operators';
import {useExtractedObservable} from '@vueuse/rxjs';

export function useUnwrappedCallMeta(
  source: () => Record<string, BehaviorSubject<any>> | undefined,
): Record<string, any> {
  const state = Vue.reactive({} as Record<string, any>);
  useExtractedObservable(source, (meta) => {
    for (const k of Object.keys(state)) delete state[k];
    const entries = Object.entries(meta).map(([name, state$]) =>
      state$.pipe(map((s) => [name, s] as const)),
    );
    return merge(...entries).pipe(
      tap(([k, val]) => {
        state[k] = val ? Vue.markRaw(val) : undefined;
      }),
    );
  });
  return state;
}
