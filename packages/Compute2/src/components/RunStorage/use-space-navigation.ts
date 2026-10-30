import * as Vue from 'vue';
import {placeName} from './run-storage';
import {SpaceIndex} from './space-index';

/** My runs, or a path of spaces from the top level; an empty path is the top level itself. */
export type RunsLocation = {myRuns: true} | {myRuns: false, path: string[]};

export const MY_RUNS: RunsLocation = {myRuns: true};
export const SPACES: RunsLocation = {myRuns: false, path: []};

export const locationKey = (location: RunsLocation) =>
  location.myRuns ? 'my' : `spaces:${location.path.join('/')}`;

export type Crumb = {label: string, depth: number};

/** Where a browser is and how it moves: into a subspace, up the way it came (also out of a linked space), or
 * to a breadcrumb. */
export function useSpaceNavigation(index: Vue.Ref<SpaceIndex | undefined>, initial: RunsLocation = MY_RUNS) {
  const location = Vue.shallowRef<RunsLocation>(initial);
  const spaceId = Vue.computed(() => location.value.myRuns ? undefined : location.value.path.at(-1));
  const space = Vue.computed(() => spaceId.value ? index.value?.nodes.get(spaceId.value)?.space : undefined);
  const top = () => index.value?.top ?? [];

  // the spaces under the location: the folders it shows, and what a search covers with everything below them
  const scope = Vue.computed(() => {
    const current = location.value;
    if (current.myRuns)
      return [];
    return spaceId.value ? [spaceId.value] : top();
  });
  const folders = Vue.computed(() => {
    const current = location.value;
    if (current.myRuns)
      return [];
    return spaceId.value ? index.value?.nodes.get(spaceId.value)?.children ?? [] : top();
  });

  const crumbs = Vue.computed<Crumb[]>(() => {
    const current = location.value;
    if (current.myRuns)
      return [{label: 'My runs', depth: 0}];
    return [{label: 'Spaces', depth: 0}, ...current.path.map((id, i) => {
      const node = index.value?.nodes.get(id);
      return {label: node ? placeName(node.space) : '…', depth: i + 1};
    })];
  });

  const go = (next: RunsLocation) => {
    if (locationKey(next) !== locationKey(location.value))
      location.value = next;
  };
  const enter = (childId: string) => {
    const current = location.value;
    if (!current.myRuns)
      go({...current, path: [...current.path, childId]});
  };
  const goTo = (depth: number) => {
    const current = location.value;
    if (!current.myRuns)
      go({...current, path: current.path.slice(0, depth)});
  };
  const up = () => {
    const current = location.value;
    if (!current.myRuns && current.path.length > 0)
      goTo(current.path.length - 1);
  };

  // a space that left the index, or stopped being under its parent, takes the path back to what is still there
  Vue.watch(index, (loaded) => {
    const current = location.value;
    if (!loaded || current.myRuns)
      return;
    const kept: string[] = [];
    for (const id of current.path) {
      const under = kept.length ? loaded.nodes.get(kept.at(-1)!)!.children : loaded.top;
      if (!loaded.nodes.has(id) || !under.includes(id))
        break;
      kept.push(id);
    }
    if (kept.length !== current.path.length)
      go({...current, path: kept});
  });

  return {location: location as Readonly<typeof location>, spaceId, space, scope, folders, crumbs, go, enter, up, goTo};
}
