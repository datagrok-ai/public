/* A stage's main panels on a twelve-column grid, laid out from a list of (panel, width): the layout is data, so a
   user's or a lab's own (the layout editor, later) is a different list, not different code. As the view narrows,
   the panels widen to half and then to the full row (container queries in studio.css), nothing hidden. */
import '../css/studio.css';
import {div} from '@datagrok-libraries/u2';

import {block} from './ui';

/** A panel's width in twelfths: a quarter, a third, a half, two thirds, the row. */
export type Span = 3 | 4 | 6 | 8 | 12;

/** A panel: what it is (title, one line on what it shows) and its content. */
export interface Panel {id: string, title: string, hint: string | HTMLElement, content: HTMLElement[]}

/** Where a panel goes: its id and its width. */
export interface Slot {id: string, span: Span}

/** The panels of `layout`, in its order and widths; a panel the layout leaves out is not shown. */
export function panelGrid(panels: Panel[], layout: Slot[], cls: string): HTMLElement {
  const byId = new Map(panels.map((p) => [p.id, p]));
  const cells = layout.filter((s) => byId.has(s.id)).map((s) => {
    const p = byId.get(s.id)!;
    const cell = div([block(p.title, p.hint, ...p.content)], 'pmx-panel');
    cell.dataset.panel = s.id;
    cell.dataset.span = String(s.span);
    return cell;
  });
  return div([div(cells, 'pmx-grid')], `pmx-grid-host ${cls}`);
}
