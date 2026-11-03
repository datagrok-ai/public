/* Moving between the app's stages: the stage's view when it is open, else a new one; Analyze with a study to analyze
   always a new one, in place of the one open. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {planView} from './plan-view';
import {analyzeView} from './analyze-view';
import {monitorView} from './monitor-view';
import type {Stage} from './ui';

const NAMES: {[stage in Stage]: string} = {plan: 'PMX Studio · Plan', monitor: 'PMX Studio · Monitor',
  analyze: 'PMX Studio · Analyze'};

export async function openStage(stage: Stage, study?: DG.DataFrame): Promise<void> {
  const open = [...grok.shell.views].find((v) => v.name === NAMES[stage]);
  // after the click that asked for it: the platform makes the clicked view current once the click is handled
  if (open && study)
    open.close();
  else if (open) {
    setTimeout(() => grok.shell.v = open);
    return;
  }
  const view = stage === 'plan' ? await planView() : stage === 'analyze' ? await analyzeView(study) :
    await monitorView();
  grok.shell.addView(view);
  // the platform puts the app's root before the path of the view it opens as the app, not of one added later
  view.path = `/apps/PmxStudio/${stage}`;
}
