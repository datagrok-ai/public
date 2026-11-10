/* The study's header, the same in Plan, Monitor and Analyze: the study's name, its groups in the colours every chart
   gives them, where the study stands (Overview · Plan — Monitor — Analyze, the stages behind it filled) and a line of
   status. The steps sit in the view's ribbon, beside the view's own actions; the study's row, over its panes. */
import '../css/studio.css';
import * as ui from 'datagrok-api/ui';
import {div, divV, span} from '@datagrok-libraries/u2';
import type {Child} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {Stage, STAGES} from './ui';
import {StudyMeta, detailsForm, emptyMeta} from './study-meta';

/** What the header says of the study: its name, its groups (the control first) and a status line; its details
 *  (ID, owner, day 0, what each group receives), edited from its name when `setMeta` keeps them. */
export interface StudyIdentity {
  name: string;
  groups: {name: string, color: string}[];
  status?: string;
  /** the study's details as they are now (read when shown and when edited: a note or an ID may have come since) */
  meta?: () => StudyMeta | null;
  setMeta?: (m: StudyMeta) => void;
  /** a group's right-click menu */
  groupMenu?: (group: string, e: MouseEvent) => void;
}

/** The header of `stage`'s view: `steps` for the ribbon (where the study stands), `row` the study over the panes;
 *  `identity` once the view knows its study, `open` brings another stage up, `owner` (the view's form) re-draws it
 *  as the study changes. */
export function studyHeader(owner: {effect: (fn: () => void) => unknown, own?: (dispose: () => void) => unknown},
  stage: Stage,
  identity: ReadonlySignal<StudyIdentity | undefined>, open: (stage: Stage) => void):
  {steps: HTMLElement, row: HTMLElement} {
  const id = span('', 'pmx-head-id');
  const name = document.createElement('button');
  name.type = 'button';
  name.className = 'pmx-head-name';
  const groups = div([], 'pmx-head-groups');
  const status = span('', 'pmx-head-status');
  const showId = (s: StudyIdentity | undefined) => id.textContent = s?.meta?.()?.id ?? '';
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
    const m = s.meta?.() ?? emptyMeta();
    const names = s.groups.length ? s.groups.map((g) => g.name) : m.groups.map((g) => g.group);
    const form = detailsForm(m, names);
    const save = ui.bigButton('Save', () => {
      const edited = form.read();
      // what the sheet does not edit stays the study's own: a note written, or an ID given, while it was open
      const now = s.meta?.();
      s.setMeta!({...edited, notes: now?.notes ?? edited.notes, id: edited.id || now?.id || ''});
      showId(s);
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
    showId(s);
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
  return {steps: div([steps(stage, open, (dispose) => owner.own?.(dispose))], 'pmx-head'), row};
}

/** A stage's page: the study's row, then the panes, which take the rest. */
export function stagePage(row: HTMLElement, panes: Child): HTMLElement {
  return divV([row, panes], 'pmx-stage');
}

/** Overview, then the three stages joined by wires; the current one marked, the ones before it filled. `own` takes
 *  what stops it watching the ribbon's size, for the view to call when it closes. */
export function steps(current: Stage, open: (stage: Stage) => void, own?: (dispose: () => void) => void): HTMLElement {
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
    const step = document.createElement(on ? 'span' : 'button');
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
  // centred in the header while the view's own buttons leave room around it, else just after them
  // (after the observer's turn: moving it within the turn would be a ResizeObserver loop)
  const place = () => {
    const ribbon = nav.closest('.d4-ribbon');
    // a view out of sight is detached: placed again when it shows
    if (!ribbon)
      return;
    fit.observe(ribbon);
    Object.assign(nav.style, {position: '', left: '', transform: ''});
    const me = nav.getBoundingClientRect();
    const own = nav.closest('.d4-ribbon-item');
    const others = [...ribbon.querySelectorAll('.d4-ribbon-item')].filter((it) => it !== own)
      .map((it) => it.getBoundingClientRect()).filter((r) => r.width > 0);
    const left = Math.max(ribbon.getBoundingClientRect().left, ...others.map((r) => r.right)) + 24;
    const box = nav.offsetParent?.getBoundingClientRect();
    if (me.left >= left || !box)
      return;
    // the view's own controls (menu, close) on the right end of the header
    const controls = ribbon.parentElement?.querySelector('.d4-view-ribbon-controls')?.getBoundingClientRect();
    const right = Math.min(box.right, controls && controls.width > 0 ? controls.left : box.right) - 8;
    // no room after the buttons either (a narrow window): in line with the header, before them
    Object.assign(nav.style, left + me.width <= right ? {left: `${left - box.left}px`, transform: 'none'} :
      {position: 'static', transform: 'none'});
  };
  const fit = new ResizeObserver(() => requestAnimationFrame(place));
  fit.observe(nav);
  own?.(() => fit.disconnect());
  return nav;
}
