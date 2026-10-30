/* confirmDialog(): asking before an action that cannot be taken back — a u2 Dialog whose footer is CANCEL and the
   action as a verb (DELETE), answering whether the user went ahead. Not `confirm`: a missing import would call the
   browser's. */
import {Dialog} from './dialog.js';
import type {DialogOkContext} from './dialog.js';
import {div, span} from '../../core/elements.js';
import {icon} from '../display/icon.js';

export interface ConfirmDialogOptions {
  title: string;
  /** The question, naming what it acts on: "Are you sure you want to delete campaign BTK follow-up (KHL-3)?" */
  message: string;
  /** A quieter line under it: what goes with it, whether it can be undone. */
  detail?: string;
  /** The action's caption, a verb: 'Delete'. 'OK' by default. */
  confirmText?: string;
  /** The action destroys something: a trash icon beside the question. */
  destructive?: boolean;
  /** The work itself, run on the action: the dialog stays open and busy until it settles, a failure stays in its
   * error row, `false` keeps it open. Without it the answer comes at the press. */
  onConfirm?: (ctx: DialogOkContext) => unknown;
  /** Blocks the page under it. */
  modal?: boolean;
  /** The automation id (`data-u2-name`). */
  name?: string;
}

/** Shows the question and answers `true` once the action ran (`onConfirm` included), `false` on CANCEL, ✕ or
 * Escape. CANCEL takes the first focus, so an Enter or a Space cancels: the platform's confirmations never focus their
 * action either (`Modal.confirm`, `Modal.confirmDelete`). */
export function confirmDialog(options: ConfirmDialogOptions): Promise<boolean> {
  return new Promise((resolve) => {
    const text = div([span(options.message, 'u2-confirm-message')], 'u2-confirm-text');
    if (options.detail !== undefined)
      text.append(span(options.detail, 'u2-confirm-detail'));
    const body = div(options.destructive === true ? [icon('trash-alt', {cls: 'u2-confirm-icon'}), text] : [text],
      'u2-confirm');
    const dialog = new Dialog(options.title, {name: options.name})
      .add(body)
      .onOK((ctx) => options.onConfirm?.(ctx), {text: options.confirmText ?? 'OK', primary: true})
      .onClose((reason) => resolve(reason === 'ok'))
      .show({modal: options.modal});
    dialog.getButton('CANCEL')?.focus();
  });
}
