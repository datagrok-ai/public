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

/** The URL path of a stage's view; the home's is the app's own. */
export const pathOf = (stage: Stage) => stage === 'home' ? '/apps/PmxStudio' : `/apps/PmxStudio/${stage}`;

export async function openStage(stage: Stage): Promise<void> {
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
    view.name = uniqueName(stage === 'home' ? 'PMX Studio' : `PMX Studio · ${STAGES.find((s) => s.stage === stage)?.label}`);
    grok.shell.addView(view);
    // the platform puts the app's root before the path of the view it opens as the app, not of one added later
    view.path = pathOf(stage);
    resolve(view);
  }));
}

/** The shell's panels a PMX Studio view hides: the context panel in every stage, and on the Overview also the
 *  browse tree, the toolbox and help, so the page stands alone. */
type Panel = 'showContextPanel' | 'showBrowse' | 'showToolbox' | 'showHelp';
const QUIET: Panel[] = ['showContextPanel'];
const QUIET_OVERVIEW: Panel[] = ['showContextPanel', 'showBrowse', 'showToolbox', 'showHelp'];

/** The state each hidden panel had before the app hid it; given back as soon as the current view no longer hides
 *  it (leaving the Overview for a stage, or leaving the app). */
const before = new Map<Panel, boolean>();
let watching = false;

/** Hides the panels above while their PMX Studio view is current (F4 still opens the context panel inside). Once
 *  per session: the app and the file handlers call it. */
export function quietContextPanel(): void {
  if (watching)
    return;
  watching = true;
  const apply = () => {
    const w = grok.shell.windows;
    const stage = stageOf(grok.shell.v);
    const hide = stage === 'home' ? QUIET_OVERVIEW : stage ? QUIET : [];
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
  grok.events.onCurrentViewChanged.subscribe(apply);
  // the app's first view may already be current when this runs
  setTimeout(apply);
}
