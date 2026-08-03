/**
 * The save / discard / cancel prompt — what stands between a user gesture and an
 * action that would silently throw away a {@link DomainFrameEditor}'s pending
 * batch.
 *
 * The editors themselves are refresh-agnostic by design: `refresh()` rebuilds and
 * discards, and deciding whether that is acceptable is the CALLER's job. This
 * module is that decision, implemented once: {@link AppView} and
 * {@link EntityListWidget} route every rebuild through {@link confirmDiscardChanges},
 * and an app that drives an editor itself should do the same.
 *
 * @module unsaved
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DomainFrameEditor} from './frame-editor';

/** What the user chose in {@link promptUnsavedChanges} (dismissing it is `'cancel'`). */
export type UnsavedOutcome = 'save' | 'discard' | 'cancel';

/** Wording of {@link promptUnsavedChanges}. */
export interface UnsavedPromptOptions {
  /** What the user is about to do — 'leave this page', 'change the filter';
   * defaults to 'continue'. */
  action?: string;
  /** Names the pending batch; defaults to the editors' table addresses. */
  subject?: string;
}

/**
 * Asks what to do with the pending changes of [editors] before something that
 * would discard them. Resolves to the user's decision WITHOUT applying it —
 * {@link confirmDiscardChanges} is the one that saves or discards.
 */
export function promptUnsavedChanges(editors: DomainFrameEditor[],
  options?: UnsavedPromptOptions): Promise<UnsavedOutcome> {
  const dirty = (editors ?? []).filter((e) => e != null && e.isDirty);
  const changes = dirty.reduce((n, e) => n + e.changeCount, 0);
  const subject = options?.subject ?? dirty.map((e) => e.table).join(', ');
  const action = options?.action ?? 'continue';
  return new Promise<UnsavedOutcome>((resolve) => {
    // Every path goes through one guard: a double-clicked button (or a close
    // notification arriving after a decision) must repeat the first answer, not
    // overwrite it.
    let decided = false;
    let closed: {unsubscribe(): void} | null = null;
    const decide = (outcome: UnsavedOutcome) => {
      // The onClose subscription outlives the dialog otherwise: its subject is
      // the platform's, not this promise's.
      closed?.unsubscribe();
      closed = null;
      if (!decided) {
        decided = true;
        resolve(outcome);
      }
    };
    const dlg = DG.Dialog.create({title: 'Unsaved changes'});
    dlg.add(ui.divV([
      ui.p(`${changes} unsaved change${changes === 1 ? '' : 's'} in ${subject}.`),
      ui.p(`Save them, discard them, or cancel and do not ${action}.`),
    ], 'ui-hint-block'));
    // Explicit indices: buttons are INSERTED into the command bar (which already
    // carries CANCEL), so the default index 0 would read DISCARD, SAVE, CANCEL.
    dlg.addButton('SAVE', () => {
      decide('save');
      dlg.close();
    }, 0);
    dlg.addButton('DISCARD', () => {
      decide('discard');
      dlg.close();
    }, 1);
    dlg.onCancel(() => decide('cancel'));
    // A dismissal (the X, Esc, a programmatic close) is a cancel — and `onClose`
    // fires for the decided paths too, where the guard keeps the first answer.
    closed = dlg.onClose.subscribe(() => decide('cancel'));
    dlg.show();
  });
}

/**
 * THE gate an app puts in front of a rebuild: resolves to whether the caller may
 * proceed, having saved or discarded whatever was pending.
 *
 * - nothing dirty — proceeds silently, no dialog;
 * - a save in flight — refuses (an editor is closed while its transaction runs,
 *   and prompting would offer a save that cannot be taken); the caller may retry
 *   after `onSavingChanged`;
 * - otherwise it prompts and applies the answer: `save` writes every pending
 *   batch (a FAILED save cancels the navigation — the user keeps their changes
 *   and the error), `discard` drops them, `cancel` stays.
 */
export async function confirmDiscardChanges(editors: DomainFrameEditor[],
  options?: UnsavedPromptOptions): Promise<boolean> {
  const all = (editors ?? []).filter((e) => e != null);
  // Mid-save first, and for EVERY editor: one with a transaction in flight
  // refuses writes, discard and refresh anyway, and prompting would offer a
  // save that cannot be taken. Retry after `onSavingChanged`.
  if (all.some((e) => e.isSaving)) {
    grok.shell.warning('Wait for the batch being saved to finish.');
    return false;
  }
  const dirty = all.filter((e) => e.isDirty);
  if (dirty.length === 0)
    return true;
  const outcome = await promptUnsavedChanges(dirty, options);
  if (outcome === 'cancel')
    return false;
  if (outcome === 'discard') {
    for (const editor of dirty)
      editor.discard();
    return true;
  }
  for (const editor of dirty)
    if (!(await editor.save()))
      return false;
  return true;
}
