/* A stage's main panels on a twelve-column grid, laid out from a list of (panel, width): the layout is data, so a
   user's or a lab's own is a different list, not different code. The layout editor (the view's menu, View > Layout…)
   orders the panels by dragging, sets their widths (a quarter, a third, a half, two thirds, the row) and hides or
   shows them; it is kept per user (grok.userSettings), the lab's default in the package settings, else the stage's
   own. As the view narrows, the panels widen to half and then to the full row (container queries in studio.css). */
import '../css/studio.css';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {div, span} from '@datagrok-libraries/u2';

import {block} from './ui';
import {_package} from './package';
import {canEditLabRules} from './outlier-rules';

/** A panel's width in twelfths: a quarter, a third, a half, two thirds, the row. */
export type Span = 3 | 4 | 6 | 8 | 12;
export const SPANS: Span[] = [3, 4, 6, 8, 12];

/** A panel: what it is (title, one line on what it shows) and its content. */
export interface Panel {id: string, title: string, hint: string | HTMLElement, content: HTMLElement[]}

/** Where a panel goes: its id and its width; `hidden`: left out. */
export interface Slot {id: string, span: Span, hidden?: boolean}

const SETTINGS = 'PmxStudio';
const keyOf = (stage: string) => `layout-${stage}`;

/** A stored layout over the stage's panels: unknown panels dropped, panels it does not name (added since) at the
 *  end as the stage places them; null: none stored, or unreadable. */
function readLayout(text: unknown, defaults: Slot[]): Slot[] | null {
  if (typeof text !== 'string' || text.trim() === '')
    return null;
  try {
    const x = JSON.parse(text);
    if (!Array.isArray(x))
      return null;
    const known = new Set(defaults.map((s) => s.id));
    const kept: Slot[] = x.filter((s) => known.has(s?.id) && SPANS.includes(s.span))
      .map((s) => ({id: s.id, span: s.span, hidden: !!s.hidden}));
    const named = new Set(kept.map((s) => s.id));
    return [...kept, ...defaults.filter((s) => !named.has(s.id))];
  } catch {
    return null;
  }
}

const labLayout = (stage: string, d: Slot[]) => readLayout(_package.settings?.[keyOf(stage)], d);
const userLayout = (stage: string, d: Slot[]) => readLayout(grok.userSettings.getValue(SETTINGS, keyOf(stage)), d);
const same = (a: Slot[], b: Slot[]) => JSON.stringify(a.map((s) => [s.id, s.span, !!s.hidden])) ===
  JSON.stringify(b.map((s) => [s.id, s.span, !!s.hidden]));

type Editor = () => void;
const editors = new WeakMap<HTMLElement, Editor>();

/** The panels of `layout`, in its order and widths (a panel it leaves out is not shown); with `stage`, the user's
 *  layout or the lab's in its place, and the grid's layout editor (editLayout). */
export function panelGrid(panels: Panel[], layout: Slot[], cls: string, stage?: string): HTMLElement {
  const cells = new Map(panels.map((p) => {
    const cell = div([block(p.title, p.hint, ...p.content)], 'pmx-panel');
    cell.dataset.panel = p.id;
    return [p.id, cell];
  }));
  const titles = new Map(panels.map((p) => [p.id, p.title]));
  const defaults = layout.filter((s) => cells.has(s.id));
  const grid = div([], 'pmx-grid');
  let current = defaults;
  const apply = (l: Slot[]) => {
    current = l;
    grid.replaceChildren(...l.filter((s) => !s.hidden).map((s) => {
      const c = cells.get(s.id)!;
      c.dataset.span = String(s.span);
      return c;
    }));
    // the charts follow their panel's new width
    window.dispatchEvent(new Event('resize'));
  };
  apply(stage ? userLayout(stage, defaults) ?? labLayout(stage, defaults) ?? defaults : defaults);
  const host = div([grid], `pmx-grid-host ${cls}`);
  if (stage)
    editors.set(host, () => layoutDialog(stage, defaults, current, titles, apply));
  return host;
}

/** Opens the layout editor of a grid panelGrid built with a stage. */
export function editLayout(grid: HTMLElement): void {
  editors.get(grid)?.();
}

/** The editor: the panels as tiles on the twelve columns, as they will sit; drag a tile to move it (or Alt+arrows on
 *  a focused one), its numbers set its width, × hides it, a hidden one comes back from the tray. */
function layoutDialog(stage: string, defaults: Slot[], now: Slot[], titles: Map<string, string>,
  apply: (l: Slot[]) => void): void {
  let work = now.map((s) => ({...s}));
  const board = div([], 'pmx-layout-board');
  const tray = div([], 'pmx-layout-tray');
  const status = span('', 'pmx-layout-status');
  let dragged: string | null = null;
  const move = (id: string, to: number) => {
    const from = work.findIndex((x) => x.id === id);
    const [m] = work.splice(from, 1);
    work.splice(Math.max(0, Math.min(work.length, to)), 0, m);
  };
  const button = (text: string, title: string, cls: string, run: () => void) => {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = cls;
    b.textContent = text;
    b.title = title;
    b.addEventListener('click', run);
    return b;
  };
  const draw = (focus?: string) => {
    const lab = labLayout(stage, defaults);
    status.textContent = same(work, lab ?? defaults) ? (lab ? 'The lab\'s layout' : 'The stage\'s own layout') :
      'Your layout';
    board.replaceChildren(...work.filter((s) => !s.hidden).map((s) => {
      const tile = div([], 'pmx-layout-tile');
      tile.style.gridColumn = `span ${s.span}`;
      tile.draggable = true;
      tile.tabIndex = 0;
      tile.dataset.panel = s.id;
      tile.setAttribute('aria-label', `${titles.get(s.id)}, ${s.span} of 12 columns`);
      const widths = div(SPANS.map((n) => button(String(n), `${n} of 12 columns`,
        `pmx-layout-w${n === s.span ? ' on' : ''}`, () => {
          s.span = n;
          draw(s.id);
        })), 'pmx-layout-widths');
      tile.append(span(titles.get(s.id) ?? s.id, 'pmx-layout-title'), widths,
        button('×', 'Hide this panel', 'pmx-layout-hide', () => {
          s.hidden = true;
          draw();
        }));
      tile.addEventListener('dragstart', (e) => {
        dragged = s.id;
        e.dataTransfer?.setData('text/plain', s.id);
        tile.classList.add('dragging');
      });
      tile.addEventListener('dragend', () => {
        dragged = null;
        tile.classList.remove('dragging');
      });
      tile.addEventListener('dragover', (e) => {
        if (dragged && dragged !== s.id)
          e.preventDefault();
      });
      tile.addEventListener('drop', (e) => {
        e.preventDefault();
        if (!dragged || dragged === s.id)
          return;
        const r = tile.getBoundingClientRect();
        const id = dragged;
        const at = work.findIndex((x) => x.id === s.id);
        const from = work.findIndex((x) => x.id === id);
        move(id, (e.clientX > r.left + r.width / 2 ? at + 1 : at) - (from < at ? 1 : 0));
        draw(id);
      });
      tile.addEventListener('keydown', (e) => {
        if (e.altKey && (e.key === 'ArrowLeft' || e.key === 'ArrowRight')) {
          e.preventDefault();
          move(s.id, work.findIndex((x) => x.id === s.id) + (e.key === 'ArrowLeft' ? -1 : 1));
          draw(s.id);
        }
      });
      return tile;
    }));
    const hidden = work.filter((s) => s.hidden);
    tray.replaceChildren(span(hidden.length ? 'Hidden:' : 'Every panel is shown', 'pmx-layout-tray-label'),
      ...hidden.map((s) => button(`+ ${titles.get(s.id) ?? s.id}`, 'Show this panel again', 'pmx-layout-chip', () => {
        s.hidden = false;
        draw(s.id);
      })));
    if (focus)
      board.querySelector<HTMLElement>(`[data-panel="${focus}"]`)?.focus();
  };
  draw();
  const help = span('Drag a panel to move it (Alt+arrows on a focused one); its numbers set its width in twelfths; ' +
    '× hides it. The header and the figures stay on top.', 'pmx-layout-help');
  const d = ui.dialog(`Layout of ${stage[0].toUpperCase()}${stage.slice(1)}`)
    .add(div([help, board, tray, status], 'pmx-layout'))
    .onOK(() => {
      // the lab's layout (or the stage's) kept as it is is no layout of one's own: a later lab change reaches it
      if (same(work, labLayout(stage, defaults) ?? defaults))
        grok.userSettings.delete(SETTINGS, keyOf(stage));
      else
        grok.userSettings.add(SETTINGS, keyOf(stage), JSON.stringify(work));
      // saved now, not at the next sync (every 10 s): a reload right after keeps it
      void grok.userSettings.flush();
      apply(work);
    });
  d.addButton('Reset', () => {
    work = (labLayout(stage, defaults) ?? defaults).map((s) => ({...s}));
    draw();
  }, undefined, 'Back to the lab\'s layout, or the stage\'s own');
  void canEditLabRules().then((can) => {
    if (!can)
      return;
    d.addButton('Save as the lab\'s', () => {
      // every setting kept: the lab's rules and group live there too
      // @ts-ignore (setSettings is declared with a Map; the platform takes the plain object, as ApiTests does)
      void _package.setSettings({..._package.settings, [keyOf(stage)]: JSON.stringify(work)}, DG.Group.allUsers)
        .then(() => {
          grok.userSettings.delete(SETTINGS, keyOf(stage));
          void grok.userSettings.flush();
          apply(work);
          grok.shell.info('Saved as the lab\'s layout');
          d.close();
        }).catch((e: any) => grok.shell.error(String(e?.message ?? e)));
    }, undefined, 'Every user\'s layout who has none of their own');
  });
  d.show({width: 760});
}
