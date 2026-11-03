import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {Observable, Subject, merge} from 'rxjs';
import {RUN_SAVED_EVENT, historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import type {RunMetadata, SaveTarget} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import type {ItemMetadata} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/view/ViewCommunication';
import {
  NESTED_ENTITY_OPTION, getRunAuthorId, getRunCreated, getRunName,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';

export const STEP_HISTORY_OPTION = 'compute2StepHistory';

/** A space, or My runs when null. */
export type Place = DG.Project | null;

export const placeName = (place: Place) => place ? (place.friendlyName || place.name) : 'My runs';

export const modelLabel = (func: DG.Func) => func.friendlyName || func.name;

export const formatSaved = (run: DG.FuncCall) => getRunCreated(run)?.format('YYYY-MM-DD HH:mm') ?? '';

const savedValue = (run: DG.FuncCall) => getRunCreated(run)?.valueOf() ?? 0;

export const bySavedDesc = (a: DG.FuncCall, b: DG.FuncCall) => savedValue(b) - savedValue(a);

const asList = (models: DG.Func | DG.Func[]) => Array.isArray(models) ? models : [models];

/** Whether a listed run belongs in the list of these models; space listings come without the package. */
export const placeRunFilter = (models: DG.Func | DG.Func[], stepHistory = false) => {
  const funcs = asList(models);
  return (run: DG.FuncCall) => funcs.some((func) => historyUtils.isRunOf(func, run)) &&
    run.options?.[NESTED_ENTITY_OPTION] !== 'true' && (!stepHistory || run.options?.[STEP_HISTORY_OPTION] === 'true');
};

/** My runs of the models, newest first. Step history keeps only runs saved as step history. */
export async function listMyModelRuns(models: DG.Func | DG.Func[], stepHistory = false): Promise<DG.FuncCall[]> {
  const runs = await historyUtils.listMyRuns(models);
  return runs.filter(placeRunFilter(models, stepHistory)).sort(bySavedDesc);
}

/** Links a run into a space under a name; a new name is saved on the run first, so every place shows it. */
export async function addToSpace(run: DG.FuncCall, space: DG.Project, name: string, replace: DG.FuncCall[] = []) {
  if (name !== getRunName(run)) {
    const edited = await historyUtils.shallowLoadRun(run.id);
    historyUtils.setRunMetadata(edited, {...historyUtils.getRunMetadata(edited), name});
    await historyUtils.updateRunMeta(edited);
  }
  await historyUtils.linkRun(run.id, space.id, replace);
}

export type SpaceAddOutcome =
  {run: DG.FuncCall, kind: 'add'} |
  {run: DG.FuncCall, kind: 'skip'} |
  {run: DG.FuncCall, kind: 'replace', replaces: DG.FuncCall[]} |
  {run: DG.FuncCall, kind: 'blocked', reason: string};

/** What adding each run under its own name does to a space; names stay unique in a space, ignoring case. */
export function planSpaceAdd(
  runs: DG.FuncCall[], listing: {items: DG.Entity[], runs: DG.FuncCall[]},
): SpaceAddOutcome[] {
  return runs.map((run): SpaceAddOutcome => {
    const name = getRunName(run);
    if (listing.runs.some((item) => item.id === run.id))
      return {run, kind: 'skip'};
    if (runs.some((other) => other.id !== run.id && historyUtils.sameRunName(getRunName(other), name)))
      return {run, kind: 'blocked', reason: `Another selected run is also named ${name}`};
    const conflict = historyUtils.nameConflictIn(listing, name);
    if (conflict.other)
      return {run, kind: 'blocked', reason: `A ${entityKind(conflict.other)} named ${name} is already there`};
    return conflict.runs.length > 0 ? {run, kind: 'replace', replaces: conflict.runs} : {run, kind: 'add'};
  });
}

/** Runs the action on each run in turn; a failure leaves the others done. */
async function forEachRun(runs: DG.FuncCall[], action: (run: DG.FuncCall, index: number) => Promise<unknown>) {
  const done: DG.FuncCall[] = [];
  const failed: {run: DG.FuncCall, error: string}[] = [];
  for (const [index, run] of runs.entries()) {
    try {
      await action(run, index);
      done.push(run);
    } catch (e: any) {
      failed.push({run, error: String(e?.message ?? e)});
    }
  }
  return {done, failed};
}

/** Links the runs a plan adds or swaps in; one failed link leaves the others in place. */
export async function addRunsToSpace(space: DG.Project, outcomes: SpaceAddOutcome[]) {
  const linked = outcomes.filter((outcome) => outcome.kind === 'add' || outcome.kind === 'replace');
  const {done, failed} = await forEachRun(linked.map((outcome) => outcome.run), (run, index) => {
    const outcome = linked[index];
    return historyUtils.linkRun(run.id, space.id, outcome.kind === 'replace' ? outcome.replaces : []);
  });
  return {added: done, failed};
}

/** Deletes runs from My runs, or removes them from a space; one failure leaves the others done. */
export const deleteOrRemoveRuns = (runs: DG.FuncCall[], space: Place) =>
  forEachRun(runs, (run) => space ? historyUtils.unlinkRun(run.id, space.id) : historyUtils.deleteRun(run));

export async function canEditSpace(place: Place): Promise<boolean> {
  return !!place && await grok.dapi.permissions.check(place, 'Edit');
}

export type RunAuthor = {id: string, name: string};

const idChunks = (ids: string[], size = 100) =>
  Array.from({length: Math.ceil(ids.length / size)}, (_, i) => ids.slice(i * size, (i + 1) * size));

const idsFilter = (ids: string[]) => `id in (${ids.map((id) => `"${id}"`).join(',')})`;

/** Authors of runs by run id. Listings come without the platform author, so runs saved before the author option
 * get it from one more request. Runs whose author isn't a known user are left out. */
export async function loadRunAuthors(runs: DG.FuncCall[]): Promise<Map<string, RunAuthor>> {
  const authorIds = new Map<string, string>();
  const unknown: string[] = [];
  for (const run of runs) {
    const id = getRunAuthorId(run);
    if (id)
      authorIds.set(run.id, id);
    else
      unknown.push(run.id);
  }
  for (const chunk of idChunks(unknown)) {
    const loaded = await grok.dapi.functions.calls.allPackageVersions()
      .filter(idsFilter(chunk)).include('session.user').list();
    for (const run of loaded) {
      const id = getRunAuthorId(run);
      if (id)
        authorIds.set(run.id, id);
    }
  }
  const names = new Map<string, string>();
  for (const chunk of idChunks([...new Set(authorIds.values())])) {
    for (const user of await grok.dapi.users.filter(idsFilter(chunk)).list())
      names.set(user.id, user.friendlyName);
  }
  const authors = new Map<string, RunAuthor>();
  for (const [runId, userId] of authorIds) {
    if (names.has(userId))
      authors.set(runId, {id: userId, name: names.get(userId)!});
  }
  return authors;
}

export const runVersion = (run: DG.FuncCall): string | undefined =>
  run.options?.['version'] || run.func?.package?.version;

/** A Save dialog result as the workflow driver takes it. */
export const asItemMetadata = ({metadata, target}: {metadata: RunMetadata, target: SaveTarget}): ItemMetadata =>
  ({title: metadata.name, description: metadata.description, tags: metadata.tags, target});

export const parseTags = (text: string) =>
  [...new Set(text.split(',').map((tag) => tag.trim()).filter((tag) => !!tag))];

/** A run's name for the Save dialog: the loaded run's name, else the workflow's title link, else model and date. */
export function defaultRunName(func: DG.Func, loadedName?: string, titleLink?: string) {
  if (loadedName)
    return loadedName;
  if (titleLink)
    return titleLink;
  return `${modelLabel(func)} ${dayjs().format('YYYY-MM-DD HH:mm')}`;
}

export function entityKind(item: DG.Entity) {
  switch (item.entityType) {
  case 'Project': return 'folder';
  case 'FuncCall': return 'run';
  case 'Script': return 'script';
  case 'DataQuery': return 'query';
  case 'DataConnection': return 'connection';
  case 'FileInfo': return 'file';
  default: return item.entityType?.toLowerCase() || 'item';
  }
}

/** What a change touched: the spaces whose content changed, and runs whose details changed wherever they are.
 * Without either, browsers reload every space. */
export type RunsChange = {spaceIds?: string[], runIds?: string[]};

const runsChangedSubject = new Subject<RunsChange>();

/** Emits after a run is saved, edited, linked or deleted, so open run browsers reload. Saves come from the platform's
 * event bus, so that a save reaches them wherever it was made. */
export const runsChanged$: Observable<RunsChange> =
  merge(runsChangedSubject, grok.events.onCustomEvent(RUN_SAVED_EVENT));

export const notifyRunsChanged = (change: RunsChange = {}) => runsChangedSubject.next(change);

export const spaceListingError = (error: string) => `Couldn't list this space: ${error}`;
