import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {FilterOption} from '../Inspector/FilterDropdown';
import {RunAuthor} from './run-storage';

export type RunFilters = {search: string, tags: string[], authors: string[]};

export const NO_FILTERS: RunFilters = {search: '', tags: [], authors: []};

/** Filter values for the caller to own, and what they do to runs; options count the runs they would keep. */
export function useRunFilters(
  runs: Vue.MaybeRefOrGetter<DG.FuncCall[]>, authors: Vue.MaybeRefOrGetter<Map<string, RunAuthor> | undefined>,
) {
  const filters = Vue.shallowRef<RunFilters>(NO_FILTERS);

  // selected tags stay listed even when the runs don't have them, so they can always be cleared
  const tagOptions = Vue.computed<FilterOption[]>(() => {
    const counts = new Map<string, number>(filters.value.tags.map((selected) => [selected, 0]));
    for (const item of Vue.toValue(runs).flatMap((run) => historyUtils.getRunMetadata(run).tags))
      counts.set(item, (counts.get(item) ?? 0) + 1);
    return [...counts].sort(([a], [b]) => a.localeCompare(b))
      .map(([value, count]) => ({value, label: value, detail: `${count}`}));
  });
  const authorOptions = Vue.computed<FilterOption[]>(() => {
    const counts = new Map<string, {name: string, count: number}>();
    const runAuthors = Vue.toValue(authors);
    for (const run of Vue.toValue(runs)) {
      const author = runAuthors?.get(run.id);
      if (author)
        counts.set(author.id, {name: author.name, count: (counts.get(author.id)?.count ?? 0) + 1});
    }
    return [...counts].sort(([, a], [, b]) => a.name.localeCompare(b.name))
      .map(([value, {name, count}]) => ({value, label: name, detail: `${count}`}));
  });

  const isActive = Vue.computed(() => {
    const {search, tags, authors: selectedAuthors} = filters.value;
    return !!search.trim() || tags.length > 0 || selectedAuthors.length > 0;
  });

  const matches = (run: DG.FuncCall) => {
    const {name, description, tags} = historyUtils.getRunMetadata(run);
    // a run must have every selected tag, as in the platform's multi-value filter
    if (filters.value.tags.some((selected) => !tags.includes(selected)))
      return false;
    // a run has one author, so selected authors are alternatives
    if (filters.value.authors.length > 0 &&
      !filters.value.authors.includes(Vue.toValue(authors)?.get(run.id)?.id ?? ''))
      return false;
    const query = filters.value.search.trim().toLowerCase();
    return !query || name.toLowerCase().includes(query) || `${description}`.toLowerCase().includes(query);
  };

  return {filters, tagOptions, authorOptions, isActive, matches};
}
