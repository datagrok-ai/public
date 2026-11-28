// Moving between the app's home and stages in one view: the stage shown replaces the one before, keeping how to
// bring it back; one with unsaved changes asks first (mayLeaveFor). Every move is an entry in the browser's history
// (the view's path, with ?plan= / ?study=), and Back and Forward walk them inside the app (route).
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import {planView} from './plan-view';
import {analyzeView} from './analyze-view';
import {monitorView} from './monitor-view';
import {hubView} from './hub-view';
import {STAGES, closeModals, plainError, stageOf} from './ui';
import {open as openSaved, tableOfRow} from './store';
import type {Stage} from './ui';
import {cruxWake} from './crux';
import {SavedPlan, parsePlan} from './plan-file';
import {dirtyViews, discardedViews, inUse, leaveAsks, letGoTables} from './native';
import {undoKey} from './panels';
import {summaryNow} from './handlers';
import {appOpened, whileAppOpen} from './app-open';
import {closeMenus} from './report';

/** Where a page of the app stands: a stage, or Compare (no stage, a page of its own). */
export type Place = Stage | 'compare';
const APP = '/apps/PmxStudio';
/** The name of every view of the app (the platform's title and the browser tab). */
const APP_NAME = 'PMX Studio';
/** The URL path of a stage's view; the home's is the app's own. */
export const pathOf = (place: Place) => place === 'home' ? APP : `${APP}/${place}`;
/** A page of the app's place: its stage, or Compare. */
export const placeOf = (v: DG.ViewBase | null | undefined): Place | undefined =>
  stageOf(v) ?? (v?.root?.dataset[STANDS_ALONE] ? 'compare' : undefined);

/** The data attribute of a view opened for an Overview row (its row id, or file path): the row brings it up again. */
export const OPENED_AS = 'pmxOpenedAs';
/** The data attribute with the saved row a stage's view shows (`plan=<id>`, `study=<id>`). */
const ROW_OF = 'pmxRowOf';

/** The view's URL path: its stage's, with the saved row it shows (?plan=<id>, ?study=<id>: what a reload, Back or a
 *  link opens, package.ts). The view the platform opened as the app gets it relative (the platform puts the app's
 *  root before it), one the app added whole. */
export function pathOfView(view: DG.ViewBase, place: Place = placeOf(view) ?? 'home'): string {
  const row = view.root.dataset[ROW_OF];
  const path = view.root.dataset[ADDED] ? pathOf(place) : place === 'home' ? '' : `/${place}`;
  return row ? `${path}?${row}` : path;
}

/** The whole address of a page of the app (what its history entry holds). */
export function addressOf(view: DG.ViewBase): string {
  const p = pathOfView(view);
  return p.startsWith(APP) ? p : APP + p;
}

/** What an address of the app stands for: its place and the saved row it names (`plan=<id>`, `study=<id>`); null
 *  for an address that is not the app's. The platform's own parameters (browse=…) are no part of it. */
export function placeAt(address: string): {place: Place, row?: string} | null {
  const url = new URL(address, 'http://x');
  const rest = url.pathname.toLowerCase();
  const app = APP.toLowerCase();
  if (rest !== app && !rest.startsWith(app + '/'))
    return null;
  const seg = rest.slice(app.length + 1).split('/')[0];
  const place = (['plan', 'monitor', 'analyze', 'compare'] as const).find((p) => p === seg) ?? 'home';
  const q = url.searchParams;
  const row = q.get('study') ? `study=${q.get('study')}` : q.get('plan') ? `plan=${q.get('plan')}` : undefined;
  return {place, row};
}

/** `address` (the platform's, with its own parameters) naming `view`'s page: its path when it names another place,
 *  and the saved row the view shows, or none (?plan= / ?study= of another page dropped). */
export function withRowOf(address: string, view: DG.ViewBase): string {
  const url = new URL(address || '/', 'http://x');
  const at = placeAt(url.pathname);
  const mine = new URL(addressOf(view), 'http://x');
  if (at?.place !== placeOf(view))
    url.pathname = mine.pathname;
  url.searchParams.delete('plan');
  url.searchParams.delete('study');
  mine.searchParams.forEach((v, k) => url.searchParams.set(k, v));
  return url.pathname + url.search + url.hash;
}

/** Whether two addresses of the app stand for the same page. */
export function samePlace(a: string, b: string): boolean {
  const x = placeAt(a);
  const y = placeAt(b);
  return !!x && !!y && x.place === y.place && x.row === y.row;
}

// ---- the history ----------------------------------------------------------------------------------------------------
// One entry per move, pushed by the app; Back and Forward replayed here, inside the app's one view (route). On Back
// the platform brings up an open view at that very path, else runs the app function again with the path (it cannot
// ask a view hosted from JS for it) and docks the view the function answers, none for null. So the app function, with
// a page of the app open, answers null and leaves the address to that page (routeOpen); the replay itself starts on
// the popstate, after the platform's (a Back between two rows of one stage runs no app function: the platform only
// brings up the view at that path). While the app moves (and while it replays), a push of the platform's (it pushes
// the current view's address on each change) rewrites the entry in hand instead of adding one. Back from the app's
// first entry leaves it as before (Datagrok's home).

/** >0 while a move is under way: the platform's pushes rewrite the entry in hand. */
let moving = 0;
/** >0 while Back or Forward is replayed: the entry exists, none is pushed. */
let replaying = 0;
/** The move being replayed (tests wait for it). */
let routing: Promise<void> = Promise.resolve();
export const routed = () => routing;
/** The last Back or Forward's address, until the platform's app function call for it (late: a server call) reads it:
 *  a call for it once the app has moved on to another page is stale and replays nothing. */
let popped: {address: string, at: number} | null = null;
/** The address the replay is to reach next (null: none waiting). */
let wanted: {address: string} | null = null;

/** Replays the moves to `address` one at a time, after the one under way: a burst of Back and Forward (each its
 *  popstate) goes to the last address only, never two replays at once (one closing the page another brings up left
 *  a page with nothing in it, QA 0.1.3 designer D3). */
function routeTo(address: string): Promise<void> {
  const idle = wanted === null;
  wanted = {address};
  if (idle) {
    routing = routing.then(async () => {
      while (wanted !== null) {
        const next = wanted;
        await route(next.address).catch((e) => console.warn('PMX Studio: a move back or forward failed', e));
        if (wanted === next)
          wanted = null;
      }
    });
  }
  return routing;
}

/** Runs `fn` as a move: the platform's pushes meanwhile, and a tick after (its view events), rewrite the entry. */
function asMove<T>(fn: () => T): T {
  moving++;
  try {
    return fn();
  } finally {
    setTimeout(() => moving--);
  }
}

/** A new entry for `address`, unless the one in hand stands for it already (the page shown brought up again). */
function pushEntry(address: string): void {
  if (samePlace(location.pathname + location.search, address))
    return;
  History.prototype.pushState.call(history, {pmx: address}, '', address);
}

/** The pushState found when the hook went in (the browser's, or another's hook): every push not turned is its. */
let wrapped: History['pushState'] | null = null;
/** Whether the hook turns pushes (a page of the app open); off, it only passes them on. */
let active = false;
/** Whether the hook turns pushes now (tests). */
export const historyHooked = (): boolean => active;
/** The hook on history.pushState, only while a page of the app is open (whileAppOpen below). It turns a push into a
 *  replace while the app moves (the platform pushes the view docked or brought up: the app pushed the entry
 *  already), or while the app's page is current and the entry in hand already stands for it (the platform may push
 *  a tick after the move). Any other push (another view current, nothing of the app moving) goes on to `wrapped`
 *  untouched. ponytail: a hook on a global, confined to the app's lifetime; a history event, when the platform
 *  gives one. */
function hook(this: History, ...args: Parameters<History['pushState']>): void {
  const own = wrapped ?? History.prototype.pushState;
  if (!active) // (the app closed while another hook sat over this one: a plain pass-through)
    return own.apply(this, args);
  const v = grok.shell.v;
  const inHand = !!v && ours(v) && samePlace(location.pathname + location.search, addressOf(v));
  // the entry already names the page shown: kept so, whatever address of the app the platform's push carries, which
  // may have lost the page's ?study= (QA 0.1.3 technician F14); an address not the app's is a view of the platform's
  // own on its way (Open as table): its entry is added, not written over the app's (F15, the URL went to /)
  if (inHand && !moving) {
    const to = String(args[2] ?? '');
    return !to || placeAt(new URL(to, location.href).pathname) ? History.prototype.replaceState.call(this, args[0],
      '', withRowOf(to || location.pathname + location.search, v!)) : own.apply(this, args);
  }
  return (moving || inHand ? History.prototype.replaceState : own).apply(this, args);
}

whileAppOpen(() => {
  // the hook goes in once and comes out only from the top: one left under another's hook (another copy of the app,
  // as the tests' bundle beside the package's, unwinding in its own order) stays there, passing through while off,
  // and is turned on again, never put over the chain that holds it (a loop) nor taken out to a null
  active = true;
  if (wrapped === null && history.pushState !== hook) {
    wrapped = history.pushState;
    history.pushState = hook;
  }
  // (with a page of the app open, current or not: a storm of Back and Forward that passed by Datagrok's home on the
  // way ends on the app's page its address names, QA 0.1.3 designer D14)
  const onPop = () => {
    const address = location.pathname + location.search;
    popped = {address, at: Date.now()};
    if (placeAt(address) && [...grok.shell.views].some((v) => ours(v)))
      void routeTo(address);
  };
  window.addEventListener('popstate', onPop);
  return () => {
    window.removeEventListener('popstate', onPop);
    active = false;
    // (another hook put over this one since: it stays, and this one only passes through)
    if (history.pushState === hook && wrapped) {
      history.pushState = wrapped;
      wrapped = null;
    }
  };
});

/** Takes the app's history in hand while a page of the app is open (quietContextPanel, showStage call it). */
export const startRouter = (): void => appOpened();

/** The app function run with a page of the app open (Back or Forward onto another stage, a link to one, the app
 *  opened again): that page takes the address and comes to the front, and true (the function answers no view: no
 *  second one). `path`: the app's path the platform passes. False with no page of the app open. */
export async function routeOpen(path?: string): Promise<boolean> {
  const open = [...grok.shell.views].find((v) => ours(v));
  if (!open)
    return false;
  // (Back or Forward: the replay began on the popstate, before the platform ran the function)
  await routing;
  // the address bar, with its row, when it is on the place asked (the platform moved it there first); else the path
  const here = location.pathname + location.search;
  const at = placeAt(here);
  const asked = path === undefined ? undefined : `${APP}/${path.replace(/^\//, '')}`;
  // the call for a Back or Forward the app has moved on from (another page shown since): stale, nothing replayed
  // (Analyze brought back over the Overview opened meanwhile). ponytail: matched by place within 10 s, as the
  // platform's call carries no token of its popstate
  const pop = popped;
  popped = null;
  if (pop && asked && Date.now() - pop.at < 10000 && placeAt(pop.address)?.place === placeAt(asked)?.place &&
    !samePlace(here, pop.address))
    return true;
  const address = at && (!asked || placeAt(asked)?.place === at.place) ? here : asked;
  if (address)
    await routeTo(address);
  else if (open.dart !== grok.shell.v?.dart) // the app opened again from elsewhere (Browse): its page as it was
    bringUp(open);
  return true;
}

/** Shows the page of the app `address` stands for, as Back or Forward reached it: the open view showing it, else the
 *  one left for another with its state, else its saved row opened, else the stage as it is. No entry is added. */
export async function route(address: string): Promise<void> {
  const at = placeAt(address);
  if (!at)
    return;
  replaying++;
  moving++;
  try {
    const cur = grok.shell.v;
    if (cur && ours(cur) && samePlace(addressOf(cur), address))
      return;
    // (a dialog of the page left goes with it)
    closeModals();
    // the page left with changes not saved asks first (Save / Discard / Cancel), rather than stay open as a second
    // page of the app behind the one Back shows (QA 0.1.3 designer D13); Cancel stays on it, its address back
    if (cur && ours(cur) && !await mayLeave([cur])) {
      pushEntry(addressOf(cur));
      return;
    }
    const open = [...grok.shell.views].find((v) => ours(v) && samePlace(addressOf(v), address));
    if (open) {
      bringUp(open);
      await new Promise((r) => setTimeout(r));
      return;
    }
    if (left.get(at.place)?.row === at.row && await returnTo(at.place))
      return;
    if (at.row) {
      const [kind, id] = at.row.split('=') as ['plan' | 'study', string];
      try {
        const v = await openRow(kind, id, at.place === 'analyze');
        await showStage(placeOf(v) as Stage ?? 'monitor', v, false);
        return;
      } catch (e) {
        grok.shell.warning(plainError(e, `The ${kind} could not be opened again.`));
      }
    }
    await openStage(at.place === 'compare' ? 'home' : at.place);
    // (a stage brought up after the click that asked for it: bringUp's tick)
    await new Promise((r) => setTimeout(r));
  } finally {
    replaying--;
    setTimeout(() => moving--);
  }
}

/** The one gate of every way out of a page: the app's views but `keep` (the page a move shows), each with changes not
 *  saved, ask first (Save / Discard / Cancel); false once one is kept. The ⌂, the steps, Follow, the demo links,
 *  another study opened and Back all come here before anything moves (QA 0.1.3 round 3: designer B1, B3; technician
 *  N1, N6, N8; analyst 7, 9; chaos 12): none leaves a "PMX Studio (2)" behind with the changes in it. */
const mayLeaveFor = (keep: DG.ViewBase) =>
  mayLeave([...grok.shell.views].filter((v) => ours(v) && v.dart !== keep.dart));

/** Each of `views` with changes not saved asks (Save / Discard / Cancel) in turn; false once one is kept (Cancel),
 *  that view then brought to the front. */
async function mayLeave(views: DG.ViewBase[]): Promise<boolean> {
  for (const v of views) {
    const ask = leaveAsks.get(v.dart);
    if (dirtyViews.get(v.dart)?.() && ask && !await ask()) {
      if (v.dart !== grok.shell.v?.dart)
        bringUp(v);
      return false;
    }
  }
  return true;
}

/** A plan or study row opened in its stage: a study with its row (so Save updates it) in Monitor, or in Analyze when
 *  finished or `analyze` (its link was Analyze's: a reload there); a plan in Plan. */
export async function openRow(kind: 'plan' | 'study', id: string, analyze = false): Promise<DG.ViewBase> {
  if (kind === 'plan')
    return planView((await openSaved('plan', id)).doc);
  const {doc, saved} = await openSaved('study', id);
  const df = tableOfRow(doc, saved, undefined, inUse);
  return analyze || doc.finished === true ? analyzeView({study: df}) : monitorView(df);
}

/** Marks `view` as showing the saved row `id`: its URL then names it (a reload opens it again, not a blank stage) and
 *  the Overview's row brings this view up rather than open another. */
export function linkRow(view: DG.ViewBase | undefined, kind: 'plan' | 'study', id: string | null | undefined): void {
  const d = view?.root?.dataset;
  const row = id ? `${kind}=${id}` : undefined;
  if (!view || !d || d[ROW_OF] === row)
    return;
  // (another table in the view, not saved yet: the URL no longer names the one before)
  if (d[ROW_OF] && d[OPENED_AS] === d[ROW_OF].split('=')[1])
    delete d[OPENED_AS];
  if (row) {
    d[ROW_OF] = row;
    d[OPENED_AS] = id!;
  } else
    delete d[ROW_OF];
  // (the entry in hand rewritten: a study saved is no move, one Back still leaves it)
  if (stageOf(view))
    asMove(() => view.path = pathOfView(view));
  syncAddress(view);
}

/** The address bar made to name `view`'s page, while it is on screen (the entry in hand rewritten, none added): the
 *  platform does not always follow a view's new path, and its own push as a view is added may carry the page
 *  before's, so a stage's URL lost its ?study= and a study not saved kept another's (QA 0.1.3 technician F14,
 *  analyst #8). */
function syncAddress(view: DG.ViewBase): void {
  const here = location.pathname + location.search;
  // (not while Back or Forward is replayed: the address bar may already be on the next entry of a burst)
  if (!replaying && wanted === null && view.dart === grok.shell.v?.dart && !samePlace(here, addressOf(view)))
    History.prototype.replaceState.call(history, {pmx: addressOf(view)}, '', addressOf(view));
}

/** The stages being opened: a second click meanwhile (a double click) waits for the same view, not another. */
const opening = new Map<string, Promise<void>>();

/** The Plan view opened with each study's plan (its text), to come back to it. */
const planViews = new Map<string, DG.ViewBase>();

/** `plan`: the plan of the study on screen (its text), when it has one: Plan opens with it, in a view of its own
 *  (the one opened with it before, while it is open); without, or one that cannot be read, Plan as it is. */
export function openStage(stage: Stage, plan?: string | null): Promise<void> {
  const key = stage === 'plan' && plan ? `plan:${plan}` : stage;
  let p = opening.get(key);
  if (!p) {
    // (the ⌂ and the steps: a dialog of the page left goes with it)
    closeModals();
    p = openStageNow(stage, stage === 'plan' ? plan : null).finally(() => opening.delete(key));
    opening.set(key, p);
  }
  return p;
}

/** What a stage view keeps when another stage replaces it: made just before it closes, it answers how to open it again
 *  as it was (null: nothing to keep). */
export type Snapshot = () => Promise<(() => Promise<DG.ViewBase>) | null>;
const snapshots = new WeakMap<ViewHost, Snapshot>();
/** Registers how the host's stage view comes back after another stage replaces it. */
export const keepOnLeaving = (host: ViewHost, snapshot: Snapshot) => void snapshots.set(host, snapshot);
/** The pages closed for another, how each opens again, the saved row it showed (its address names it: Back finds
 *  it by it) and (Plan) the study's plan it was opened with. */
type Left = {reopen: () => Promise<DG.ViewBase>, row?: string, plan?: string};
const left = new Map<Place, Left>();
/** Whether `stage` has something to come back to: its view open, or closed for another with its state kept. */
export const canReturn = (stage: Stage) => left.has(stage) || [...grok.shell.views].some((v) => stageOf(v) === stage);

/** Closes every view of the app but `keep`, each stage keeping how to come back (keepOnLeaving); one with changes
 *  not saved stays open. */
async function leaveOthers(keep: DG.ViewBase): Promise<void> {
  for (const v of [...grok.shell.views]) {
    if (v.dart === keep.dart || !ours(v) || dirtyViews.get(v.dart)?.())
      continue;
    const place = placeOf(v);
    const snapshot = place && place !== 'home' ? snapshots.get(ViewHost.of(v)) : undefined;
    let reopen: Awaited<ReturnType<Snapshot>> = null;
    // (changes discarded: not kept to come back to, Back opens the study as saved; QA 0.1.3 round 3 chaos 1)
    if (place && snapshot && !discardedViews.has(v.dart)) {
      try {
        reopen = await snapshot();
      } catch (e) {
        console.warn('PMX Studio: a stage left without its state', e);
      }
    }
    // (brought to the front again meanwhile, by a later move: it stays, and is not kept as left either; closed, the
    // platform showed it again with nothing in it, QA 0.1.3 designer D3)
    if (v.dart === grok.shell.v?.dart)
      continue;
    if (place && reopen) {
      left.set(place, {reopen, row: v.root.dataset[ROW_OF],
        plan: [...planViews].find(([, pv]) => pv.dart === v.dart)?.[0]});
    }
    // (gone meanwhile: a second move while the first was keeping its state)
    if ([...grok.shell.views].some((x) => x.dart === v.dart))
      v.close();
  }
  // (its plain name back once the one it replaced is gone)
  if ([...grok.shell.views].some((x) => x.dart === keep.dart))
    keep.name = uniqueName(APP_NAME, keep);
  // (a view closed lets the current object go: the study of the one kept is put back)
  setTimeout(summaryNow);
}

/** An open view of the app brought to the front in place of the others, after the click that asked for it: the
 *  platform makes the clicked view current once the click is handled, which would undo a switch made in it. */
export function bringUp(view: DG.ViewBase): void {
  // (Back or Forward bringing it up: its entry exists)
  const push = !replaying;
  setTimeout(async () => {
    const gone = () => ![...grok.shell.views].some((v) => v.dart === view.dart);
    // (closed meanwhile, by a move made since: made current, the platform would dock it again with nothing in it,
    // QA 0.1.3 designer D3; the page the address names is shown instead)
    if (gone()) {
      void routeTo(location.pathname + location.search);
      return;
    }
    // a page with changes not saved asks before it goes; Cancel: it stays, the address as it was (chaos 13)
    if (!await mayLeaveFor(view) || gone())
      return;
    asMove(() => {
      if (push)
        pushEntry(addressOf(view));
      grok.shell.v = view;
      void leaveOthers(view);
    });
    setTimeout(() => syncAddress(view));
  });
}

/** The view shown for `stage` from what it kept when another stage replaced it; what it kept is let go only once
 *  it opens again (one that fails says so and stays, to try again). */
async function returnTo(stage: Place): Promise<boolean> {
  const l = left.get(stage);
  if (!l)
    return false;
  let reopened: DG.ViewBase;
  try {
    reopened = await l.reopen();
  } catch (e) {
    grok.shell.error(plainError(e, `${STAGES.find((s) => s.stage === stage)?.label ?? 'Compare'} could not open again. ` +
      'Try again.'));
    return true;
  }
  // (its address names its study from the first: a step to it lost ?study= until its first reading, QA 0.1.3
  // technician F14)
  if (l.row) {
    const [kind, id] = l.row.split('=') as ['plan' | 'study', string];
    linkRow(reopened, kind, id);
  }
  const view = await showStage(stage === 'compare' ? 'home' : stage, reopened, false);
  // (a page kept by its Cancel: what this one kept stays, to come back to)
  if (view.dart !== reopened.dart)
    return true;
  left.delete(stage);
  if (l.plan)
    planViews.set(l.plan, view);
  return true;
}

async function openStageNow(stage: Stage, plan?: string | null): Promise<void> {
  if (plan) {
    // by the platform's view, not the JS wrapper: grok.shell.views may hand back another wrapper of the same view
    const kept = planViews.get(plan);
    const own = kept && [...grok.shell.views].find((v) => v.dart === kept.dart);
    if (own) {
      bringUp(own);
      return;
    }
    if (left.get('plan')?.plan === plan && await returnTo('plan'))
      return;
    let saved: SavedPlan | null = null;
    try {
      saved = parsePlan(plan);
    } catch {} // a damaged plan: Analyze and Monitor say so; Plan opens as it is
    if (saved) {
      planViews.set(plan, await showStage('plan', await planView(saved), false));
      return;
    }
  }
  const open = [...grok.shell.views].find((v) => stageOf(v) === stage);
  // after the click that asked for it: the platform makes the clicked view current once the click is handled
  if (open) {
    bringUp(open);
    return;
  }
  if (await returnTo(stage))
    return;
  // a stage's first opening takes a moment (its engine, its data): said in the task bar meanwhile
  const progress = DG.TaskBarProgressIndicator.create(`Opening ${STAGES.find((s) => s.stage === stage)?.label}...`);
  try {
    const view = stage === 'plan' ? await planView() : stage === 'analyze' ? await analyzeView() :
      stage === 'monitor' ? await monitorView() : await hubView();
    // (a move between the stages of the study in hand, or to the Overview and back: what the others left is kept)
    await showStage(stage, view, false);
  } finally {
    progress.close();
  }
}

/** `base`, else `base (2)`, … : the shell keys its list of views by name (a second "PMX Studio" logs "Mapping for PMX
 *  Studio exists" each time the list is drawn). */
export function uniqueName(base: string, self?: DG.ViewBase): string {
  const taken = new Set([...grok.shell.views].filter((v) => v.dart !== self?.dart).map((v) => v.name));
  let name = base;
  for (let k = 2; taken.has(name); k++)
    name = `${base} (${k})`;
  return name;
}

/** A stage's view opened with something in it (a study, a plan): always a view of its own, in place of the app's
 *  other views (leaveOthers); `fresh`: another study than the one in hand, so what the stages left behind is
 *  forgotten (a plan followed into Monitor is still kept: it leaves after this). Added after the click that asked
 *  for it: the platform makes the clicked view current once the click is handled (in a task of its own), so a view
 *  ready at once (its data at hand) would be put behind the one clicked. */
export async function showStage(stage: Stage, view: DG.ViewBase, fresh = true, keepUnsaved = false):
  Promise<DG.ViewBase> {
  // any page in place of one with changes not saved (another study, or a stage of the same: the ⌂, a step, Follow):
  // it asks first (Save / Discard / Cancel), rather than stay open behind as "PMX Studio (2)" (QA 0.1.3 technician
  // F13, F16; round 3 designer B1); Cancel keeps it, and the new page is let go. `keepUnsaved`: a page showing the
  // very table not saved (Monitor's Analyze this study), whose changes go on with it: the one left stays open
  if (!keepUnsaved && !await mayLeaveFor(view)) {
    // (and its tables with it, at once: Cancel changes nothing, no "X (2)" left in the Study field)
    letGoTables(view);
    ViewHost.of(view).dispose();
    const others = [...grok.shell.views].filter((v) => ours(v) && v.dart !== view.dart);
    return others.find((v) => dirtyViews.get(v.dart)?.()) ?? others[0];
  }
  if (fresh)
    left.clear();
  appOpened();
  // (Back or Forward showing it: its entry exists)
  const push = !replaying;
  return new Promise((resolve) => setTimeout(() => asMove(() => {
    // every page of the app is "PMX Studio" (the study is named in its header): "(2)" only while the page it
    // replaces closes (leaveOthers), or beside one kept for its unsaved changes
    view.name = uniqueName(APP_NAME, view);
    // (added views come pinned, whatever the user did: the panels read the app's pin instead, see applyPanels)
    view.root.dataset[ADDED] = 'true';
    const place = placeOf(view) ?? stage;
    // the entry first: a push of the platform's as the view is added then rewrites it, not the one before
    if (push)
      pushEntry(addressOf(view));
    grok.shell.addView(view);
    // the platform puts the app's root before the path of the view it opens as the app, not of one added later
    view.path = pathOfView(view, place);
    void leaveOthers(view);
    resolve(view);
  })));
}

/** The shell's panels a PMX Studio view hides: the context panel in every view of the app (it opens on a click on a
 *  group or an animal: ui.ts showDepth), and the Browse tree only while the view is pinned, for the page to have the
 *  whole screen. Opening the app, or coming back to the Overview, leaves Browse, the toolbox and help as the user had
 *  them (a one-click preview of the app used to close Browse). */
export type Panel = 'showContextPanel' | 'showBrowse';
/** The data attribute (camel-cased) of a page of the app that is no stage (Compare): treated as the Overview. */
export const STANDS_ALONE = 'pmxAlone';

/** The panels hidden while a view is current: `stage` its stage (null: none), `alone` a page of the app that is no
 *  stage, `pinned` the view pinned. A study's stage keeps the context panel: its study is there (handlers.ts). */
export function quietPanels(stage: string | null, alone: boolean, pinned: boolean): Panel[] {
  if (!stage && !alone)
    return [];
  const study = !!stage && stage !== 'home';
  return [...study ? [] : ['showContextPanel' as const], ...pinned ? ['showBrowse' as const] : []];
}

const ours = (v: DG.ViewBase | null) => !!stageOf(v) || !!v?.root?.dataset[STANDS_ALONE];
/** The data attribute of a view the app added itself (showStage): the platform adds every such view pinned. */
export const ADDED = 'pmxAdded';

/** Shows and hides the panels (`w`: grok.shell.windows) for the current view `v`: the ones it hides are hidden, their
 *  state kept in `before`; the ones hidden before and no longer hidden (unpinned, another PMX view, leaving the app)
 *  given back as they were. Pinned is the user's: a view the platform opened (the app from Browse, a file) is pinned
 *  when the user opened it for good (a double click) or pinned its preview, and that is kept in `app`; a view the
 *  app added (a stage, Compare) is always pinned by the platform, so it follows `app`. Whether `v` counts as pinned. */
export function applyPanels(w: {[p in Panel]: boolean}, v: DG.ViewBase | null, before: Map<Panel, boolean>,
  app: {pinned: boolean}): boolean {
  if (ours(v) && !v!.root.dataset[ADDED])
    app.pinned = !!v!.isPinned;
  const pinned = ours(v) && app.pinned;
  const hide = quietPanels(stageOf(v) ?? null, !!v?.root?.dataset[STANDS_ALONE], pinned);
  for (const [p, on] of [...before]) {
    if (!hide.includes(p)) {
      w[p] = on;
      before.delete(p);
    }
  }
  for (const p of hide) {
    if (!before.has(p))
      before.set(p, w[p]);
    w[p] = false;
  }
  return pinned;
}

/** The state each hidden panel had before the app hid it; given back as soon as the current view no longer hides
 *  it (unpinned, or leaving the app). */
const before = new Map<Panel, boolean>();

/** Hides the panels above while their PMX Studio view is current (F4 still opens the context panel inside), with the
 *  app's other listeners on the session: all only while a page of the app is open (app-open.ts). The app and the
 *  file handlers call it. */
export const quietContextPanel = (): void => appOpened();

whileAppOpen(() => {
  const computes = (v: DG.ViewBase | null) => !!stageOf(v) && stageOf(v) !== 'home'; // (Plan, Monitor, Analyze)
  // the app's pin (see applyPanels), and the current view's isPinned as last read: a pin changes no view, so it is
  // watched (below)
  const app = {pinned: false};
  let pinnedNow = false;
  const apply = () => {
    if (computes(grok.shell.v)) // a stage opened or come back to: the engine at rest wakes
      cruxWake();
    pinnedNow = !!grok.shell.v?.isPinned;
    applyPanels(grok.shell.windows, grok.shell.v, before, app);
  };
  // and again a tick later: back on a pinned Overview from a stage, the platform shows its Browse tree after the event
  const viewSub = grok.events.onCurrentViewChanged.subscribe(() => {
    apply();
    setTimeout(apply, 0);
    // (a menu of the page left does not stay open over the page shown: File opened by itself as a page came up, QA
    // 0.1.3 designer D17)
    // (and the context panel shows the study of the page brought up, not the one left: QA 0.1.3 round 3 designer
    // prev 18)
    if (ours(grok.shell.v)) {
      setTimeout(closeMenus);
      setTimeout(summaryNow);
    }
  });
  // the platform has no event for a view pinned (DG.View.fromRoot views get no afterPersist call, onViewChanged
  // does not fire on a pin): its state is read once a second while a page of the app is current, and the panels
  // change only when it changes (the user may open Browse again on a pinned view)
  // ponytail: a 1 s poll of one getter, only while the app is open; an event, when the platform gives one
  const poll = setInterval(() => {
    if (!document.hidden && ours(grok.shell.v) && !!grok.shell.v?.isPinned !== pinnedNow)
      apply();
  }, 1000);
  // the tab back on a stage that computes: the engine at rest wakes
  const onVisible = () => {
    if (!document.hidden && computes(grok.shell.v))
      cruxWake();
  };
  // a field of a stage that computes taken in hand: the value follows a second or more later, the sets open meanwhile
  const onFocus = (e: FocusEvent) => {
    if (computes(grok.shell.v) && grok.shell.v!.root.contains(e.target as Node))
      cruxWake();
  };
  // Ctrl+Z on a page of the app, outside a field and a grid: the platform's own undo has nothing of ours to take
  // back, and runs the app again (a second Overview; QA 0.1.3 technician #11). Kept from it; a page being customized
  // takes it first (panels.ts, the same capture on the document)
  const onKey = (e: KeyboardEvent) => {
    // Esc closes a menu left open over a page of the app (the ribbon's File stayed open after Ctrl+O or a day's
    // right click, and took no Esc: QA 0.1.3 designer D17)
    if (e.key === 'Escape' && ours(grok.shell.v) && closeMenus()) {
      e.preventDefault();
      e.stopPropagation();
      return;
    }
    if (ours(grok.shell.v) && undoKey(e, true) && !(e.target as HTMLElement | null)?.closest?.('.d4-grid')) {
      e.preventDefault();
      e.stopPropagation();
    }
  };
  // a click elsewhere closes a ribbon menu left open over a page of the app (designer prev 17)
  const onPointer = (e: PointerEvent) => {
    const at = e.target as HTMLElement | null;
    if (e.isTrusted && ours(grok.shell.v) && !at?.closest?.('[role="menubar"], .d4-menu-popup'))
      closeMenus();
  };
  document.addEventListener('visibilitychange', onVisible);
  document.addEventListener('focusin', onFocus);
  document.addEventListener('keydown', onKey, true);
  document.addEventListener('pointerdown', onPointer, true);
  // the app's first view may already be current when this runs
  const first = setTimeout(apply);
  return () => {
    viewSub.unsubscribe();
    clearInterval(poll);
    clearTimeout(first);
    document.removeEventListener('visibilitychange', onVisible);
    document.removeEventListener('focusin', onFocus);
    document.removeEventListener('keydown', onKey, true);
    document.removeEventListener('pointerdown', onPointer, true);
    // (the panels the app hid given back: no page of the app is current now)
    applyPanels(grok.shell.windows, null, before, app);
  };
});
