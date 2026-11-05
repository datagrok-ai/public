import * as Vue from 'vue';
import {EMPTY, Observable, Subject, defer, of} from 'rxjs';
import {catchError, concatMap} from 'rxjs/operators';
import {useSubscription} from '@vueuse/rxjs';
import {RunsChange, runsChanged$} from './run-storage';
import {SpaceIndex, loadSpaceIndex, reloadSpaces, spacesWithRun} from './space-index';

type IndexUpdate = (loaded?: SpaceIndex) => Observable<SpaceIndex | undefined>;

/** The spaces of a tree, loaded when it mounts and kept up to date as runs change. */
export function useSpaceIndex() {
  const index = Vue.shallowRef<SpaceIndex>();
  const error = Vue.ref<string>();
  const updates = new Subject<IndexUpdate>();

  // one update at a time, each starting from the index the previous ones left; unmounting stops the one running
  useSubscription(updates.pipe(concatMap((next) => defer(() => next(index.value)).pipe(catchError((e) => {
    error.value = String(e?.message ?? e);
    return EMPTY;
  })))).subscribe((loaded) => {
    if (loaded) {
      index.value = loaded;
      error.value = undefined;
    }
  }));
  const load = () => updates.next(() => loadSpaceIndex());

  // only the spaces a change touched are listed again; an unknown change reloads them all
  const touched = (loaded: SpaceIndex, change: RunsChange) => {
    if (!change.spaceIds && !change.runIds)
      return undefined;
    const withRuns = (change.runIds ?? []).flatMap((id) => spacesWithRun(loaded, id));
    return [...new Set([...change.spaceIds ?? [], ...withRuns])];
  };

  Vue.onMounted(load);
  useSubscription(runsChanged$.subscribe((change) => updates.next((loaded) => {
    const spaceIds = loaded ? touched(loaded, change) : undefined;
    if (!loaded || !spaceIds)
      return loadSpaceIndex();
    return spaceIds.length > 0 ? reloadSpaces(loaded, spaceIds) : of(undefined);
  })));

  return {index: index as Readonly<typeof index>, error: error as Readonly<typeof error>, reload: load};
}
