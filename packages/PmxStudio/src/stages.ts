/* Moving between the app's home and stages: the stage's view when it is open, else a new one. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {planView} from './plan-view';
import {analyzeView} from './analyze-view';
import {hubView} from './hub-view';
import type {Stage} from './ui';

export const NAMES: {[stage in Stage]?: string} = {home: 'PMX Studio', plan: 'PMX Studio · Plan',
  analyze: 'PMX Studio · Analyze'};

/** The URL path of a stage's view; the home's is the app's own. */
export const pathOf = (stage: Stage) => stage === 'home' ? '/apps/PmxStudio' : `/apps/PmxStudio/${stage}`;

export async function openStage(stage: Stage): Promise<void> {
  const open = [...grok.shell.views].find((v) => v.name === NAMES[stage]);
  // after the click that asked for it: the platform makes the clicked view current once the click is handled
  if (open) {
    setTimeout(() => grok.shell.v = open);
    return;
  }
  const view = stage === 'plan' ? await planView() : stage === 'analyze' ? await analyzeView() :
    stage === 'home' ? await hubView() : null;
  if (view)
    showStage(stage, view);
}

/** A stage's view opened with something in it (a study, a plan): always a view of its own. */
export function showStage(stage: Stage, view: DG.ViewBase): void {
  grok.shell.addView(view);
  // the platform puts the app's root before the path of the view it opens as the app, not of one added later
  view.path = pathOf(stage);
}
