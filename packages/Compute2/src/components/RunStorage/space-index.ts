import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {EMPTY, Observable, defer, from, of} from 'rxjs';
import {catchError, expand, filter, map, mergeMap, toArray} from 'rxjs/operators';

/** A space with its subspaces and the runs linked into it; `error` says why it couldn't be listed. */
export type SpaceNode = {space: DG.Project, children: string[], runs: DG.FuncCall[], error?: string};

/** The spaces list and every space below it, with its subspaces and the runs linked into it. A subspace linked
 * into another space is the same node under both parents, and `parents` lists them. A subspace shared without its
 * parent isn't reachable from the spaces list, so it isn't here. */
export type SpaceIndex = {top: string[], nodes: Map<string, SpaceNode>, parents: Map<string, string[]>};

export type RunFilter = (run: DG.FuncCall) => boolean;

export type SpaceListing = {subspaces: DG.Project[], runs: DG.FuncCall[]};
export type ListSpace = (space: DG.Project) => Observable<SpaceListing>;

const MAX_CONCURRENT_REQUESTS = 4;

const nameOf = (space: DG.Project) => (space.friendlyName || space.name).toLowerCase();

export const sortSpaces = (spaces: DG.Project[]) => {
  const unique = new Map(spaces.filter((space) => space.isSpace).map((space) => [space.id, space]));
  return [...unique.values()].sort((a, b) => nameOf(a).localeCompare(nameOf(b)));
};

// one request lists the owned and linked subspaces and the linked runs
export const listSpace: ListSpace = (space) =>
  defer(() => grok.dapi.spaces.id(space.id).children.ofTypes('Project,FuncCall', true).list()).pipe(map((items) => ({
    subspaces: sortSpaces(items.filter((item): item is DG.Project => item instanceof DG.Project)),
    runs: items.filter((item): item is DG.FuncCall => item.entityType === 'FuncCall'),
  })));

type WalkStep = {space: DG.Project} | {node: SpaceNode};

/** Lists the spaces and everything below them, a few at a time and each once: a space linked under several parents,
 * or one `skip` keeps, isn't listed again. A space that can't be listed is a node with the error, and the walk goes
 * on. Unsubscribing stops it. */
export function walkSpaces(spaces: DG.Project[], list: ListSpace = listSpace,
  skip: (id: string) => boolean = () => false, concurrency = MAX_CONCURRENT_REQUESTS): Observable<SpaceNode> {
  const queued = new Set<string>();
  const steps = (candidates: DG.Project[]): WalkStep[] => candidates
    .filter((space) => !skip(space.id) && !queued.has(space.id) && !!queued.add(space.id))
    .map((space) => ({space}));
  // the index is kept in shallow refs, and the entities' Dart handles must not be proxied
  const visit = (space: DG.Project) => list(space).pipe(
    mergeMap(({subspaces, runs}): WalkStep[] => [{node: {space: Vue.markRaw(space),
      children: subspaces.map((sub) => sub.id), runs: runs.map((run) => Vue.markRaw(run))}}, ...steps(subspaces)]),
    catchError((e) => of<WalkStep>({node: {space: Vue.markRaw(space), children: [], runs: [],
      error: String(e?.message ?? e)}})));
  return from(steps(spaces)).pipe(
    expand((step) => 'space' in step ? visit(step.space) : EMPTY, concurrency),
    filter((step): step is {node: SpaceNode} => 'node' in step),
    map((step) => step.node));
}

export function makeIndex(top: string[], nodes: Map<string, SpaceNode>): SpaceIndex {
  const parents = new Map<string, string[]>();
  for (const [id, node] of nodes) {
    for (const child of node.children)
      parents.set(child, [...parents.get(child) ?? [], id]);
  }
  return {top, nodes, parents};
}

export const loadSpaceIndex = (list: ListSpace = listSpace): Observable<SpaceIndex> =>
  defer(() => grok.dapi.spaces.list()).pipe(mergeMap((spaces) => {
    const top = sortSpaces(spaces);
    return walkSpaces(top, list).pipe(toArray(), map((nodes) =>
      makeIndex(top.map((space) => space.id), new Map(nodes.map((node) => [node.space.id, node])))));
  }));

/** Lists these spaces again, and any new subspaces they now have; spaces nothing links to any more stay until the
 * next full load. A space that can't be listed again keeps what it had, with the error. */
export function reloadSpaces(index: SpaceIndex, spaceIds: string[],
  list: ListSpace = listSpace): Observable<SpaceIndex> {
  const relisted = new Set(spaceIds.filter((id) => index.nodes.has(id)));
  const spaces = [...relisted].map((id) => index.nodes.get(id)!.space);
  return walkSpaces(spaces, list, (id) => index.nodes.has(id) && !relisted.has(id)).pipe(toArray(), map((listed) => {
    const nodes = new Map(index.nodes);
    for (const node of listed) {
      const old = nodes.get(node.space.id);
      nodes.set(node.space.id, node.error && old ? {...old, error: node.error} : node);
    }
    return makeIndex(index.top, nodes);
  }));
}

export type LocatedRun = {run: DG.FuncCall, spaceId: string};

/** The matching runs of the spaces and of everything below them, each run once, with the space it is in. */
export function subtreeRuns(index: SpaceIndex, spaceIds: string[], filter: RunFilter): LocatedRun[] {
  const visited = new Set<string>();
  const seen = new Set<string>();
  const found: LocatedRun[] = [];
  const visit = (id: string) => {
    const node = index.nodes.get(id);
    if (!node || visited.has(id))
      return;
    visited.add(id);
    for (const run of node.runs) {
      if (!seen.has(run.id) && filter(run)) {
        seen.add(run.id);
        found.push({run, spaceId: id});
      }
    }
    node.children.forEach(visit);
  };
  spaceIds.forEach(visit);
  return found;
}

/** The number of matching runs in each space and below it. */
export const subtreeCounts = (index: SpaceIndex, filter: RunFilter) =>
  new Map([...index.nodes.keys()].map((id) => [id, subtreeRuns(index, [id], filter).length]));

/** The spaces a space is a subspace of. */
export const parentsOf = (index: SpaceIndex, spaceId: string) =>
  (index.parents.get(spaceId) ?? []).flatMap((id) => index.nodes.get(id)?.space ?? []);

/** The favorite spaces of an index and every space above them. */
export function favoriteSpaces(index: SpaceIndex, favorites: Set<string>): Set<string> {
  const kept = new Set<string>();
  const keep = (id: string) => {
    if (kept.has(id))
      return;
    kept.add(id);
    index.parents.get(id)?.forEach(keep);
  };
  [...favorites].filter((id) => index.nodes.has(id)).forEach(keep);
  return kept;
}

/** The spaces a run is linked into. */
export const spacesWithRun = (index: SpaceIndex, runId: string) =>
  [...index.nodes.entries()].filter(([, node]) => node.runs.some((run) => run.id === runId)).map(([id]) => id);
