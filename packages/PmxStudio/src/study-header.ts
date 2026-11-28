/* The study's header, the same in Plan, Monitor and Analyze: the study's name, its groups in the colours every chart
   gives them, where the study stands (Overview · Plan — Monitor — Analyze, the stages behind it filled) and a line of
   status. The steps sit in the view's ribbon, beside the view's own actions; the study's row, over its panes. */
import '../css/studio.css';
import * as ui from 'datagrok-api/ui';
import {div, divV, span} from '@datagrok-libraries/u2';
import type {Child} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {okNamed, showCentred, Stage, STAGES, modal, objectPanel, swatch} from './ui';
import {StudyMeta, detailsForm, emptyMeta} from './study-meta';
import {myEntries} from './store';
import {PmxDay, bindObject, showSummary, studySummaries} from './handlers';
import type {PmxGroup} from './handlers';

/** The groups the header shows as dots (a control and 15 arms, the most a plan has); the others counted. */
const HEAD_DOTS = 16;

/** What the header says of the study: its name, its groups (the control first) and a status line; its details
 *  (ID, owner, day 0, what each group receives), edited from its name when `setMeta` keeps them. */
export interface StudyIdentity {
  name: string;
  groups: {name: string, color: string}[];
  status?: string;
  /** more on the status line, on hover */
  statusTip?: () => string;
  /** the study's details as they are now (read when shown and when edited: a note or an ID may have come since) */
  meta?: () => StudyMeta | null;
  setMeta?: (m: StudyMeta) => void;
  /** a blinded study before its groups are revealed: the details sheet does not show what each group receives */
  groupsHidden?: boolean;
  /** a study not saved yet (read when shown): its ID is only the one it will ask for (its plan's, its file's), which
   *  the store may refuse and replace, so the header says "Not saved yet" and names that ID on hover */
  unsaved?: () => boolean;
  /** a group as the platform's object (handlers.ts): its right click offers its actions */
  groupObject?: (group: string, color: string) => PmxGroup;
}

/** The header of `stage`'s view: `steps` for the ribbon (PMX Studio, and where the study stands), `row`
 *  the study over the panes; `identity` once the view knows its study, `open` brings another stage up (the Overview
 *  too), `owner` (the view's form) re-draws it as the study changes; `details` opens the study's details, as its name
 *  does. `off`: the stages that have nothing to show now, each with the reason (greyed, said on hover). */
export function studyHeader(owner: {effect: (fn: () => void) => unknown, own?: (dispose: () => void) => unknown},
  stage: Stage,
  identity: ReadonlySignal<StudyIdentity | undefined>, open: (stage: Stage) => void,
  off?: ReadonlySignal<StepsOff>):
  {steps: HTMLElement, row: HTMLElement, details: (then?: () => void) => void, refresh: () => void} {
  const id = span('', 'pmx-head-id');
  const name = document.createElement('button');
  name.type = 'button';
  name.className = 'pmx-head-name';
  const groups = div([], 'pmx-head-groups');
  const status = span('', 'pmx-head-status');
  // (a copy says whose the original is, as the Overview does)
  const showId = (s: StudyIdentity | undefined) => {
    const m = s?.meta?.();
    const unsaved = !!s?.unsaved?.();
    id.textContent = [unsaved ? 'Not saved yet' : m?.id ?? '', m?.copyOf ? `copy of ${m.copyOf}` : '']
      .filter((x) => x).join(' · ');
    id.title = unsaved ? `The study ID is given when it is saved${m?.id ? `: ${m.id} if no other study has it, ` +
      'else a new one' : ''}.` : '';
  };
  // the details in a dialog wide enough for the groups' schedule grid; `then` once it closes, saved or not
  const details = (then?: () => void) => {
    const s = identity.peek();
    if (!s?.setMeta)
      return;
    const m = s.meta?.() ?? emptyMeta();
    const names = s.groups.length ? s.groups.map((g) => g.name) : [...new Set(m.groups.map((g) => g.group))];
    // the IDs the user's plans and studies have: another's is refused (the list as the Overview reads it)
    const taken = myEntries().then((es) => es.map((e) => e.id));
    const form = detailsForm(m, names, s.groupsHidden, taken);
    const dialog = modal('study-details', 'Study details');
    if (!dialog)
      return;
    dialog.add(form.root).onOK(() => {
      if (form.problem() !== null)
        return;
      const edited = form.read();
      // what the dialog does not edit stays the study's own: a note written, or an ID given, while it was open
      const now = s.meta?.();
      s.setMeta!({...edited, notes: now?.notes ?? edited.notes, id: edited.id || now?.id || '',
        colors: now?.colors ?? edited.colors});
      showId(s);
    });
    const save = okNamed(dialog, 'Save');
    const allow = () => save.disabled = form.problem() !== null;
    form.onProblem(allow);
    allow();
    if (then)
      dialog.onClose.subscribe(() => then());
    showCentred(dialog, {width: Math.min(900, window.innerWidth - 40)});
  };
  name.addEventListener('click', () => details());
  owner.effect(() => {
    const s = identity.value;
    name.textContent = s?.name ?? '';
    name.disabled = !s?.setMeta;
    name.title = s?.setMeta ? 'The study\'s details: ID, owner, day 0, what each group receives' : '';
    showId(s);
    // past HEAD_DOTS groups the rest as "+N more", their names on hover: 50 dots overflowed the header and cut the
    // summary (QA 0.1.3 A13)
    const all = s?.groups ?? [];
    const rest = all.slice(HEAD_DOTS);
    const more = rest.length ? span(`+${rest.length} more`, 'pmx-head-more') : null;
    if (more) {
      Object.assign(more.style, {fontSize: '12px', color: 'var(--dg-text-color-muted, #7d8792)', whiteSpace: 'nowrap'});
      more.title = rest.map((g) => g.name).join('\n');
    }
    groups.replaceChildren(...all.slice(0, HEAD_DOTS).map((g) => {
      // a dot in the group's colour, its name on hover: the study's name keeps its room
      const chip = span('', 'pmx-head-group');
      chip.style.background = g.color;
      chip.title = g.name;
      chip.setAttribute('aria-label', g.name);
      if (s?.groupObject)
        bindObject(s.groupObject(g.name, g.color), chip);
      return chip;
    }), ...more ? [more] : []);
    status.textContent = s?.status ?? '';
    status.title = s?.statusTip?.() ?? '';
  });
  const row = div([div([name, id, groups, status], 'pmx-head-study')], 'pmx-head-block');
  owner.effect(() => row.hidden = !identity.value);
  // the study in the context panel while nothing else is picked there (handlers.ts), drawn again as it changes
  studySummaries.set(row, () => {
    const s = identity.peek();
    return !s ? null : new PmxDay(s.name, () => objectPanel({title: s.name,
      // (as the header: no ID before the study is saved, which gives it one)
      subtitle: [s.unsaved?.() ? 'Not saved yet' : s.meta?.()?.id, s.status].filter((x) => x).join(' · '),
      // (a blinded study's groups stay hidden here too)
      facts: s.groupsHidden ? [] : [{label: 'Groups', value: divV(s.groups.map((g) => swatch(g.name, g.color)))}]}));
  });
  owner.effect(() => {
    void identity.value;
    showSummary(row);
  });
  // the study named at the bar's right end (its details on a click, as its name over the panes)
  const chip = document.createElement('button');
  chip.type = 'button';
  chip.className = 'pmx-head-chip';
  chip.addEventListener('click', () => details());
  const showChip = (s: StudyIdentity | undefined) => {
    chip.textContent = s ? [s.name, s.unsaved?.() ? '' : s.meta?.()?.id ?? ''].filter((x) => x).join(' · ') : '';
    chip.hidden = !s;
    chip.disabled = !s?.setMeta;
    chip.title = name.title;
  };
  owner.effect(() => showChip(identity.value));
  return {steps: headBar(div([crumb(() => open('home')), steps(stage, open, off, (fn) => owner.effect(fn))],
    'pmx-head-center'), chip), row, details, refresh: () => {
    showId(identity.peek());
    showChip(identity.peek());
    // (the context panel too: after Save as it kept the old name and "Not saved yet", QA 0.1.3 round 3 chaos 15)
    showSummary(row);
  }};
}

/** The bar in the view's ribbon: the ⌂ and the steps centred, `right` (the study, or what Compare compares) at its
 *  right end; the ribbon's menu stays on its left. */
export function headBar(center: HTMLElement, right: HTMLElement): HTMLElement {
  const bar = div([center, right], 'pmx-head');
  // centred on the whole top bar (the view's ribbon row), not on the room the menu leaves it: the grid centres the
  // steps in the bar, then they move left by the difference, never past the bar's start (over the menu)
  let row: Element | null = null;
  // the ribbon's menu (File, View, Study, Help): the platform may lay it over the start of the bar, so the bar keeps
  // clear of it (its left padding, set only when it changes: a change re-runs this) and the steps never go under it,
  // at any width or study name; narrower, the study is cut first (QA 0.1.3 round 3: designer B2, B6; technician N17)
  // (the menu is found again whenever it is not on the page, and each one found is watched: the platform puts it in
  // after the bar, or anew, with no change of the bar's size; only the first was watched, and the ⌂ and the steps
  // stayed under the menu laid in later, QA 0.1.3 round 4 designer B2)
  let menu: Element | null = null;
  const watched = new WeakSet<Element>();
  // (each watched for its size and for its items put in, renamed or folded with no change of any size: the menu's own
  // box may stay as it was while its items run past it, and the steps sat under Study and Help at 1600 px, QA 0.1.3
  // browser check)
  const watch = (e: Element | null) => {
    if (e && !watched.has(e)) {
      watched.add(e);
      seen.observe(e);
      changes.observe(e, {childList: true, subtree: true, characterData: true});
    }
  };
  const place = () => {
    row ??= bar.closest('.d4-ribbon')?.parentElement?.closest('.d4-ribbon') ?? null;
    if (!row)
      return;
    if (!menu?.isConnected)
      menu = row.querySelector('[role="menubar"]') ?? row.parentElement?.querySelector('[role="menubar"]') ?? null;
    watch(row);
    watch(menu);
    const [r, b] = [row, bar].map((e) => e.getBoundingClientRect());
    const over = Math.ceil(ribbonEnd(row, menu, bar) + 12 - b.left);
    const pad = over > 0 && over < b.width ? `${over}px` : '';
    if (bar.style.paddingLeft !== pad)
      bar.style.paddingLeft = pad;
    center.style.translate = '';
    const c = center.getBoundingClientRect();
    const start = b.left + Math.max(0, over);
    const shift = Math.min(0, Math.max((r.left + r.right - c.left - c.right) / 2, start - c.left));
    center.style.translate = shift ? `${Math.round(shift)}px` : '';
  };
  // (on the next frame: observing and padding inside the callback is a "ResizeObserver loop" error)
  let frame = 0;
  const soon = () => {
    if (!bar.isConnected && row)
      return void window.removeEventListener('resize', soon);
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(place);
  };
  const seen = new ResizeObserver(soon);
  const changes = new MutationObserver(soon);
  seen.observe(bar);
  seen.observe(center);
  // (and the window resized)
  window.addEventListener('resize', soon);
  return bar;
}

/** Where the ribbon's own items end left of `bar` (px from the window's left; -Infinity: none there): the furthest
 *  right of the visible items (labels, icons: elements holding no other) of the ribbon and its menu in the bar's line
 *  and left of its middle, the bar and what holds it left out. Each item measured, not the menu's box, which its
 *  items may run past (QA 0.1.3 browser check: elementFromPoint(560, 16) was the ribbon's "Study"). */
export function ribbonEnd(row: Element, menu: Element | null, bar: HTMLElement): number {
  const b = bar.getBoundingClientRect();
  let end = -Infinity;
  for (const scope of menu && !row.contains(menu) ? [row, menu] : [row]) {
    for (const e of scope.querySelectorAll('*')) {
      // (an empty box, a spacer, is no item)
      if (e.childElementCount > 0 || bar.contains(e) || !(e.textContent ?? '').trim() &&
        !(e instanceof SVGElement || e.tagName === 'I' || e.tagName === 'IMG'))
        continue;
      const x = e.getBoundingClientRect();
      if (x.width > 0 && x.height > 0 && x.bottom > b.top && x.top < b.bottom && x.left < b.left + b.width / 2 &&
        getComputedStyle(e).visibility !== 'hidden')
        end = Math.max(end, x.right);
    }
  }
  return end;
}

/** The stages a step leads to that have nothing to show now, each with why (null or absent: it is on). */
export type StepsOff = {[stage: string]: string | null};

/** The home icon before the steps, leading to the Overview; `here`, after it, on a page that is no stage (Compare). */
export function crumb(home: () => void, here?: string): HTMLElement {
  const back = document.createElement('button');
  back.type = 'button';
  back.className = 'pmx-crumb-home';
  back.append(ui.iconFA('home'));
  back.title = 'Overview';
  back.setAttribute('aria-label', 'Overview');
  back.addEventListener('click', home);
  return div([back, ...here ? [span(here, 'pmx-crumb-here')] : []], 'pmx-crumb');
}

/** A stage's page: the study's row, then the panes, which take the rest. */
export function stagePage(row: HTMLElement, panes: Child): HTMLElement {
  return divV([row, panes], 'pmx-stage');
}

/** The three stages joined by wires; the current one marked, the ones before it filled. A stage in `off` (followed
 *  by `follow`, the view's effect) is greyed and takes no click, its reason on hover. */
export function steps(current: Stage, open: (stage: Stage) => void,
  off?: ReadonlySignal<StepsOff>, follow?: (fn: () => void) => unknown): HTMLElement {
  const nav = div([], 'pmx-head-steps');
  nav.setAttribute('role', 'navigation');
  nav.setAttribute('aria-label', 'Study stages');
  const list = STAGES.filter((s) => s.stage !== 'home');
  const at = list.findIndex((s) => s.stage === current);
  const buttons = new Map<string, {el: HTMLElement, tip: string}>();
  list.forEach(({stage, label, tip}, i) => {
    if (i > 0)
      nav.append(span('', `pmx-head-wire${i <= at ? ' lit' : ''}`));
    const on = stage === current;
    const step = document.createElement(on ? 'span' : 'button');
    step.className = `pmx-head-step${on ? ' on' : i < at ? ' past' : ''}`;
    step.dataset.stage = stage;
    // (on a narrow window, 1280 px and the like, a step other than the current one is its dot, named on hover:
    // studio.css) ponytail: by the window's width, not the ribbon's room; a long name with the context panel open
    // may still crowd it
    step.append(span('', 'pmx-head-node'), span(label, 'pmx-head-label'));
    step.setAttribute('aria-label', label);
    step.title = tip;
    if (on)
      step.setAttribute('aria-current', 'step');
    if (step instanceof HTMLButtonElement) {
      step.type = 'button';
      step.addEventListener('click', () => {
        if (!step.classList.contains('is-off'))
          open(stage);
      });
      buttons.set(stage, {el: step, tip});
    }
    nav.append(step);
  });
  const greyed = (why: StepsOff) => buttons.forEach(({el, tip}, stage) => {
    const reason = why[stage] ?? null;
    el.classList.toggle('is-off', !!reason);
    el.setAttribute('aria-disabled', String(!!reason));
    // (its name first: on a narrow screen the step is its dot alone, studio.css)
    el.title = `${el.getAttribute('aria-label')}: ${reason ?? tip}`;
  });
  if (off && follow)
    follow(() => greyed(off.value));
  else
    greyed({});
  return nav;
}
