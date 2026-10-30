import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {EMPTY, from} from 'rxjs';
import {watchExtractedObservable} from '@vueuse/rxjs';
import {RunAuthor, loadRunAuthors} from './run-storage';

/** The authors of the runs, by run id, loaded again when the runs change; the last load wins. */
export function useRunAuthors(runs: Vue.MaybeRefOrGetter<DG.FuncCall[]>) {
  const authors = Vue.shallowRef<Map<string, RunAuthor>>();
  const error = Vue.ref<string>();

  watchExtractedObservable(() => Vue.toValue(runs), (items) => {
    authors.value = undefined;
    error.value = undefined;
    return items.length > 0 ? from(loadRunAuthors(items)) : EMPTY;
  }, (loaded) => authors.value = loaded, {
    onError: (e: any) => error.value = `Couldn't load the authors: ${e?.message ?? e}`,
  }, {immediate: true});

  return {authors: authors as Readonly<typeof authors>, error: error as Readonly<typeof error>};
}
