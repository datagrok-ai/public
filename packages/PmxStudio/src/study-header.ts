/* The study's header, the same in Plan, Monitor and Analyze: the study's name, its groups in the colours every chart
   gives them, where the study stands (Overview · Plan — Monitor — Analyze, the stages behind it filled) and a line of
   status. It takes the place of the stage strip at the top of the view; the view's actions stay in its ribbon. */
import '../css/studio.css';
import {div, span} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {MONITOR_READY, Stage, STAGES} from './ui';
import {StudyMeta, detailsDialog, emptyMeta} from './study-meta';

/** What the header says of the study: its name, its groups (the control first) and a status line; its details
 *  (ID, owner, day 0, what each group receives), edited from its name when `setMeta` keeps them. */
export interface StudyIdentity {
  name: string;
  groups: {name: string, color: string}[];
  status?: string;
  meta?: StudyMeta | null;
  setMeta?: (m: StudyMeta) => void;
}

/** The header of `stage`'s view: `identity` once the view knows its study, `open` brings another stage up, `owner`
 *  (the view's form) re-draws it as the study changes. */
export function studyHeader(owner: {effect: (fn: () => void) => unknown}, stage: Stage,
  identity: ReadonlySignal<StudyIdentity | undefined>, open: (stage: Stage) => void): HTMLElement {
  const id = span('', 'pmx-head-id');
  const name = document.createElement('button');
  name.type = 'button';
  name.className = 'pmx-head-name';
  const groups = div([], 'pmx-head-groups');
  const status = span('', 'pmx-head-status');
  // the details as last saved here: the identity's own are read again only when the study changes
  let meta: StudyMeta | null = null;
  const showId = () => id.textContent = meta?.id ?? '';
  name.addEventListener('click', () => {
    const s = identity.peek();
    if (!s?.setMeta)
      return;
    const m = meta ?? emptyMeta();
    const names = s.groups.length ? s.groups.map((g) => g.name) : m.groups.map((g) => g.group);
    detailsDialog(s.name, m, names, (saved) => {
      meta = saved;
      showId();
      s.setMeta!(saved);
    });
  });
  owner.effect(() => {
    const s = identity.value;
    name.textContent = s?.name ?? '';
    name.disabled = !s?.setMeta;
    name.title = s?.setMeta ? 'The study\'s details: ID, owner, day 0, what each group receives' : '';
    meta = s?.meta ?? null;
    showId();
    groups.replaceChildren(...(s?.groups ?? []).map((g) => {
      const chip = span('', 'pmx-head-group');
      const box = span('', 'pmx-head-swatch');
      box.style.background = g.color;
      chip.append(box, span(g.name));
      return chip;
    }));
    status.textContent = s?.status ?? '';
  });
  const study = div([id, name, groups, status], 'pmx-head-study');
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
