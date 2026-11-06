/* The study's header, the same in Plan, Monitor and Analyze: the study's name, its groups in the colours every chart
   gives them, where the study stands (Overview · Plan — Monitor — Analyze, the stages behind it filled) and a line of
   status. It takes the place of the stage strip at the top of the view; the view's actions stay in its ribbon. */
import '../css/studio.css';
import * as ui from 'datagrok-api/ui';
import {div, divV, span} from '@datagrok-libraries/u2';
import type {Child} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {MONITOR_READY, Stage, STAGES} from './ui';
import {StudyMeta, detailsForm, emptyMeta} from './study-meta';

/** What the header says of the study: its name, its groups (the control first) and a status line; its details
 *  (ID, owner, day 0, what each group receives), edited from its name when `setMeta` keeps them. */
export interface StudyIdentity {
  name: string;
  groups: {name: string, color: string}[];
  status?: string;
  meta?: StudyMeta | null;
  setMeta?: (m: StudyMeta) => void;
  /** a group's right-click menu */
  groupMenu?: (group: string, e: MouseEvent) => void;
}

/** The header of `stage`'s view: `steps` for the ribbon (where the study stands), `row` the study over the panes;
 *  `identity` once the view knows its study, `open` brings another stage up, `owner` (the view's form) re-draws it
 *  as the study changes. */
export function studyHeader(owner: {effect: (fn: () => void) => unknown}, stage: Stage,
  identity: ReadonlySignal<StudyIdentity | undefined>, open: (stage: Stage) => void):
  {steps: HTMLElement, row: HTMLElement} {
  const id = span('', 'pmx-head-id');
  const name = document.createElement('button');
  name.type = 'button';
  name.className = 'pmx-head-name';
  const groups = div([], 'pmx-head-groups');
  const status = span('', 'pmx-head-status');
  // the details as last saved here: the identity's own are read again only when the study changes
  let meta: StudyMeta | null = null;
  const showId = () => id.textContent = meta?.id ?? '';
  // the details open under the header, the study still in view: a sheet, not a dialog; again, or Escape, closes it
  const sheet = div([], 'pmx-head-sheet');
  sheet.hidden = true;
  const close = () => {
    sheet.hidden = true;
    sheet.replaceChildren();
    name.setAttribute('aria-expanded', 'false');
  };
  sheet.addEventListener('keydown', (e) => {
    if (e.key === 'Escape')
      close();
  });
  name.addEventListener('click', () => {
    const s = identity.peek();
    if (!s?.setMeta || !sheet.hidden)
      return close();
    const m = meta ?? emptyMeta();
    const names = s.groups.length ? s.groups.map((g) => g.name) : m.groups.map((g) => g.group);
    const form = detailsForm(m, names);
    const save = ui.bigButton('Save', () => {
      meta = form.read();
      showId();
      s.setMeta!(meta);
      close();
    });
    sheet.replaceChildren(form.root, div([save, ui.button('Cancel', close)], 'pmx-head-sheet-actions'));
    sheet.hidden = false;
    name.setAttribute('aria-expanded', 'true');
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
      if (s?.groupMenu) {
        chip.addEventListener('contextmenu', (e) => {
          // the ribbon has a menu of its own
          e.preventDefault();
          e.stopPropagation();
          s.groupMenu!(g.name, e);
        });
      }
      return chip;
    }));
    status.textContent = s?.status ?? '';
  });
  const row = div([div([name, id, groups, status], 'pmx-head-study'), sheet], 'pmx-head-block');
  owner.effect(() => row.hidden = !identity.value);
  return {steps: div([steps(stage, open)], 'pmx-head'), row};
}

/** A stage's page: the study's row, then the panes, which take the rest. */
export function stagePage(row: HTMLElement, panes: Child): HTMLElement {
  return divV([row, panes], 'pmx-stage');
}

/** Overview, then the three stages joined by wires; the current one marked, the ones before it filled. */
export function steps(current: Stage, open: (stage: Stage) => void): HTMLElement {
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
