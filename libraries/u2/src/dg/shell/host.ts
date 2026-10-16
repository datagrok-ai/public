/* A u2 control placed into platform chrome — a docked panel, an accordion pane, a widget the
   context panel shows — goes when the platform kills that element, as a view's content does. */
import * as DG from 'datagrok-api/dg';
import type {Control} from '../../core/component.js';

/** `control.root`, wired to the platform's kill channel: closing the docked panel, the view around
 * it or the pane it sits in disposes `control` and everything it owns. Chrome the platform never
 * kills — an element kept across views — stays the builder's to dispose. */
export function host(control: Control): HTMLElement {
  const root = control.root;
  // core contract (2026-08): a docked pane's ✕ kills elements carrying this attribute
  root.setAttribute('data-kill-on-close', 'true');
  DG.Widget.registerCleanup(root, () => control.dispose());
  return root;
}
