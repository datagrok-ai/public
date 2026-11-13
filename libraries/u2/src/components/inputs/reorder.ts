/* Rows a user puts in order — RecordsInput's and ChecklistInput's: a drag by the row's handle, or Alt+↑/↓ from
   anywhere in the row. The handle and the row being dragged wear css/reorder.css. */
import {on} from '../../core/elements.js';
import {icon} from '../display/icon.js';

export class Reorder {
  /** A row's drag handle (`u2-reorder-grip`, part `handle`). */
  static handle(): HTMLElement {
    const handle = icon('grip-vertical', {cls: 'u2-reorder-grip'});
    handle.dataset.u2Part = 'handle';
    return handle;
  }

  /** Wires `row` — built in its own row scope, so the listeners go with it. `rows()` answers the rows in their
   * order now; `move` puts one elsewhere (live while a drag runs), `done` follows the gesture's last move. */
  static wire(row: HTMLElement, handle: HTMLElement, rows: () => HTMLElement[],
    move: (from: number, to: number) => void, done: () => void): void {
    const shift = (to: (from: number, count: number) => number): boolean => {
      const list = rows();
      const from = list.indexOf(row);
      const target = Math.max(0, Math.min(list.length - 1, to(from, list.length)));
      if (from < 0 || target === from)
        return false;
      move(from, target);
      return true;
    };
    on<KeyboardEvent>(row, 'keydown', (e) => {
      if (!e.altKey || e.defaultPrevented || (e.key !== 'ArrowUp' && e.key !== 'ArrowDown'))
        return;
      e.preventDefault();
      const active = document.activeElement as HTMLElement | null;
      if (shift((from) => from + (e.key === 'ArrowUp' ? -1 : 1)))
        done();
      // a moved element loses the focus it held; the move lands when the handler's batch ends
      queueMicrotask(() => {
        if (active !== null && active !== document.activeElement && active.isConnected)
          active.focus();
      });
    });
    let dragging = false;
    on<PointerEvent>(handle, 'pointerdown', (e) => {
      if (e.button !== 0)
        return;
      e.preventDefault();
      handle.setPointerCapture?.(e.pointerId);
      dragging = true;
      row.classList.add('u2-reorder-dragging');
    });
    on<PointerEvent>(handle, 'pointermove', (e) => {
      if (dragging)
        shift(() => rows().filter((el) => el !== row && Reorder._middle(el) < e.clientY).length);
    });
    const end = () => {
      if (!dragging)
        return;
      dragging = false;
      row.classList.remove('u2-reorder-dragging');
      done();
    };
    on(handle, 'pointerup', end);
    on(handle, 'pointercancel', end);
  }

  private static _middle(el: HTMLElement): number {
    const rect = el.getBoundingClientRect();
    return rect.top + rect.height / 2;
  }
}
