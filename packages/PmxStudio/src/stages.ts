/* Moving between the app's home and stages: the stage's view when it is open, else a new one; a stage opened with
   something in it (a study, a plan) is a view of its own. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {planView} from './plan-view';
import {analyzeView} from './analyze-view';
import {monitorView} from './monitor-view';
import {hubView} from './hub-view';
import {STAGES, stageOf} from './ui';
import type {Stage} from './ui';
import {cruxWake} from './crux';
import {SavedPlan, parsePlan} from './plan-file';

/** The URL path of a stage's view; the home's is the app's own. */
export const pathOf = (stage: Stage) => stage === 'home' ? '/apps/PmxStudio' : `/apps/PmxStudio/${stage}`;

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
    p = openStageNow(stage, stage === 'plan' ? plan : null).finally(() => opening.delete(key));
    opening.set(key, p);
  }
  return p;
}

async function openStageNow(stage: Stage, plan?: string | null): Promise<void> {
  if (plan) {
    // by the platform's view, not the JS wrapper: grok.shell.views may hand back another wrapper of the same view
    const kept = planViews.get(plan);
    const own = kept && [...grok.shell.views].find((v) => v.dart === kept.dart);
    if (own) {
      setTimeout(() => grok.shell.v = own);
      return;
    }
    let saved: SavedPlan | null = null;
    try {
      saved = parsePlan(plan);
    } catch {} // a damaged plan: Analyze and Monitor say so; Plan opens as it is
    if (saved) {
      planViews.set(plan, await showStage('plan', await planView(saved)));
      return;
    }
  }
  const open = [...grok.shell.views].find((v) => stageOf(v) === stage);
  // after the click that asked for it: the platform makes the clicked view current once the click is handled
  if (open) {
    setTimeout(() => grok.shell.v = open);
    return;
  }
  // a stage's first opening takes a moment (its engine, its data): said in the task bar meanwhile
  const progress = DG.TaskBarProgressIndicator.create(`Opening ${STAGES.find((s) => s.stage === stage)?.label}...`);
  try {
    const view = stage === 'plan' ? await planView() : stage === 'analyze' ? await analyzeView() :
      stage === 'monitor' ? await monitorView() : await hubView();
    await showStage(stage, view);
  } finally {
    progress.close();
  }
}

/** `base`, else `base (2)`, … : the shell keys its list of views by name (a second "PMX Studio" logs "Mapping for PMX
 *  Studio exists" each time the list is drawn). */
export function uniqueName(base: string): string {
  const taken = new Set([...grok.shell.views].map((v) => v.name));
  let name = base;
  for (let k = 2; taken.has(name); k++)
    name = `${base} (${k})`;
  return name;
}

/** A stage's view opened with something in it (a study, a plan): always a view of its own. Added after the click
 *  that asked for it: the platform makes the clicked view current once the click is handled (in a task of its own),
 *  so a view ready at once (its data at hand) would be put behind the one clicked. */
export function showStage(stage: Stage, view: DG.ViewBase): Promise<DG.ViewBase> {
  return new Promise((resolve) => setTimeout(() => {
    const label = STAGES.find((s) => s.stage === stage)?.label;
    view.name = uniqueName(stage === 'home' ? 'PMX Studio' : `PMX Studio · ${label}`);
    grok.shell.addView(view);
    // the platform puts the app's root before the path of the view it opens as the app, not of one added later
    view.path = pathOf(stage);
    resolve(view);
  }));
}

/** The shell's panels a PMX Studio view hides: the context panel in every view of the app (it opens on a click on a
 *  group or an animal: ui.ts showDepth), and the Browse tree only while the view is pinned, for the page to have the
 *  whole screen. Opening the app, or coming back to the Overview, leaves Browse, the toolbox and help as the user had
 *  them (a one-click preview of the app used to close Browse). */
type Panel = 'showContextPanel' | 'showBrowse';
/** The data attribute (camel-cased) of a page of the app that is no stage (Compare): treated as the Overview. */
export const STANDS_ALONE = 'pmxAlone';

/** The panels hidden while a view is current: `stage` its stage (null: none), `alone` a page of the app that is no
 *  stage, `pinned` the view pinned. */
export function quietPanels(stage: string | null, alone: boolean, pinned: boolean): Panel[] {
  if (!stage && !alone)
    return [];
  return pinned ? ['showContextPanel', 'showBrowse'] : ['showContextPanel'];
}

/** The state each hidden panel had before the app hid it; given back as soon as the current view no longer hides
 *  it (unpinned, or leaving the app). */
const before = new Map<Panel, boolean>();
let watching = false;

/** Hides the panels above while their PMX Studio view is current (F4 still opens the context panel inside). Once
 *  per session: the app and the file handlers call it. */
export function quietContextPanel(): void {
  if (watching)
    return;
  watching = true;
  const computes = (v: DG.ViewBase | null) => !!stageOf(v) && stageOf(v) !== 'home'; // (Plan, Monitor, Analyze)
  const ours = (v: DG.ViewBase | null) => !!stageOf(v) || !!v?.root?.dataset[STANDS_ALONE];
  // the current view's pinned state as last applied: a pin changes no view, so it is watched (below)
  let pinnedNow = false;
  const apply = () => {
    const w = grok.shell.windows;
    const v = grok.shell.v;
    if (computes(v)) // a stage opened or come back to: the engine at rest wakes
      cruxWake();
    pinnedNow = ours(v) && !!v?.isPinned;
    const hide = quietPanels(stageOf(v) ?? null, !!v?.root?.dataset[STANDS_ALONE], pinnedNow);
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
  };
  // and again a tick later: back on a pinned Overview from a stage, the platform shows its Browse tree after the event
  grok.events.onCurrentViewChanged.subscribe(() => {
    apply();
    setTimeout(apply, 0);
  });
  // the platform has no event for a view pinned (DG.View.fromRoot views get no afterPersist call): its state is read
  // once a second, and the panels change only when it changes (the user may open Browse again on a pinned view)
  // ponytail: a 1 s poll of one getter; an event, when the platform gives one
  setInterval(() => {
    const v = grok.shell.v;
    if ((ours(v) && !!v?.isPinned) !== pinnedNow)
      apply();
  }, 1000);
  // the tab back on a stage that computes: the engine at rest wakes
  document.addEventListener('visibilitychange', () => {
    if (!document.hidden && computes(grok.shell.v))
      cruxWake();
  });
  // a field of a stage that computes taken in hand: the value follows a second or more later, the sets open meanwhile
  document.addEventListener('focusin', (e) => {
    if (computes(grok.shell.v) && grok.shell.v!.root.contains(e.target as Node))
      cruxWake();
  });
  // the app's first view may already be current when this runs
  setTimeout(apply);
}
