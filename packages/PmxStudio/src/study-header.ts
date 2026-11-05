/* The study's header, the same in Plan, Monitor and Analyze: the study's name, its groups in the colours every chart
   gives them, where the study stands (Overview · Plan — Monitor — Analyze, the stages behind it filled) and a line of
   status. It takes the place of the stage strip at the top of the view; the view's actions stay in its ribbon. */
import '../css/studio.css';
import {div, span} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {MONITOR_READY, Stage, STAGES} from './ui';

/** What the header says of the study: its name, its groups (the control first) and a status line. */
export interface StudyIdentity {
  name: string;
  groups: {name: string, color: string}[];
  status?: string;
}

/** The header of `stage`'s view: `identity` once the view knows its study, `open` brings another stage up, `owner`
 *  (the view's form) re-draws it as the study changes. */
export function studyHeader(owner: {effect: (fn: () => void) => unknown}, stage: Stage,
  identity: ReadonlySignal<StudyIdentity | undefined>, open: (stage: Stage) => void): HTMLElement {
  const name = span('', 'pmx-head-name');
  const groups = div([], 'pmx-head-groups');
  const status = span('', 'pmx-head-status');
  owner.effect(() => {
    const s = identity.value;
    name.textContent = s?.name ?? '';
    groups.replaceChildren(...(s?.groups ?? []).map((g) => {
      const chip = span('', 'pmx-head-group');
      const box = span('', 'pmx-head-swatch');
      box.style.background = g.color;
      chip.append(box, span(g.name));
      return chip;
    }));
    status.textContent = s?.status ?? '';
  });
  const study = div([name, groups, status], 'pmx-head-study');
  return div([steps(stage, open), study], 'pmx-head');
}

/** Overview, then the three stages joined by wires; the current one marked, the ones before it filled. */
function steps(current: Stage, open: (stage: Stage) => void): HTMLElement {
  const nav = div([], 'pmx-head-steps');
  nav.setAttribute('role', 'navigation');
  nav.setAttribute('aria-label', 'Study stages');
  const at = STAGES.findIndex((s) => s.stage === current);
  STAGES.forEach(({stage, label, tip}, i) => {
    if (i === 1)
      nav.append(span('·', 'pmx-head-dot'));
    else if (i > 1)
      nav.append(span('', `pmx-head-wire${i <= at ? ' lit' : ''}`));
    const on = stage === current;
    const ready = stage !== 'monitor' || MONITOR_READY;
    const step = document.createElement(on || !ready ? 'span' : 'button');
    step.className = `pmx-head-step${on ? ' on' : i > 0 && i < at ? ' past' : ''}`;
    if (stage !== 'home')
      step.append(span('', 'pmx-head-node'));
    step.append(span(label));
    step.title = tip;
    if (on)
      step.setAttribute('aria-current', 'step');
    if (step instanceof HTMLButtonElement) {
      step.type = 'button';
      step.addEventListener('click', () => open(stage));
    }
    nav.append(step);
  });
  return nav;
}
