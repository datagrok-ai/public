/* Save, Discard and New over a session — a source's own session of one today, the phase-2
   session over several tables tomorrow. State ⇒ disabled: nothing to save, or a save in flight;
   Save also while the session's validity blocks it;
   New is the one permission-gated button (hidden without `insert`), the fields and actions
   having degraded by access already. Controls rather than bare buttons, so the effects they run
   have an owner — `appView` disposes ribbon controls with the content. */
import {Control, button} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';
import type {DomainSession} from '../core/session.js';
import type {RowView} from '@datagrok-libraries/u2';

/** A session, or a source standing for its own. */
export type SessionTarget = DomainSession | DomainSource;

export interface SessionButtonOptions {
  text?: string;
}

export class SessionButton extends Control {
  readonly button: HTMLButtonElement;
  readonly session: DomainSession;

  constructor(target: SessionTarget, text: string,
    options: {primary?: boolean, enabled?: () => boolean, run: (session: DomainSession) => unknown}) {
    const session = SessionButton.sessionOf(target);
    const el = button(text, () => void options.run(session), {primary: options.primary});
    super(el);
    this.button = el;
    this.session = session;
    this.effect(() => el.disabled = session.isSaving.value || !(options.enabled?.() ?? session.isDirty.value));
  }

  static sessionOf(target: SessionTarget): DomainSession {
    return target instanceof DomainSource ? target.session : target;
  }
}

/** Save: the batch as one transaction through the session, which announces what it saved and
 * hands the focus back to the paired form; the session's primary button, where Tab from a form's
 * last field lands. */
export function saveButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const session = SessionButton.sessionOf(target);
  const control = new SessionButton(target, options.text ?? 'Save', {primary: true,
    enabled: () => session.isDirty.value && session.validity.value === null, run: () => session.save()});
  control.root.dataset.u2 = 'save-button';
  control.button.title = 'Save (Ctrl+S)';
  session.primaryButton = control.button;
  control.own(() => {
    if (session.primaryButton === control.button)
      session.primaryButton = undefined;
  });
  return control;
}

/** `leave`: while `when` holds, Discard is on over a clean session and runs `leave.run` instead — off an
 * untouched create page. */
export function discardButton(target: SessionTarget,
  options: SessionButtonOptions & {leave?: {when: () => boolean, run: () => void}} = {}): SessionButton {
  const leave = options.leave;
  const control = new SessionButton(target, options.text ?? 'Discard', {
    enabled: leave === undefined ? undefined : () => SessionButton.sessionOf(target).isDirty.value || leave.when(),
    run: (session) => leave?.when() ? leave.run() : session.discard()});
  control.root.dataset.u2 = 'discard-button';
  control.button.title = 'Discard changes';
  return control;
}

/** "New": a draft over `values` made current — for a source, the row a create form or a list then
 * edits; a function receives the row that was current (the one just saved, say) so an app carries
 * fields over from one entry to the next. The press is a change of its own — the draft counts
 * against the session from the start, unlike the create page's initial row. Permission ⇒ hidden:
 * absent without the `insert` capability; state ⇒ disabled until the table is loaded. */
export function newButton(source: DomainSource,
  values: Record<string, unknown> | ((last: RowView | null) => Record<string, unknown>) = {},
  options: SessionButtonOptions = {}): Control {
  const el = button(options.text ?? 'New', () => source.newRow(
    typeof values === 'function' ? values(source.currentRow.peek()) : values));
  const control = new Control(el);
  control.root.dataset.u2 = 'new-button';
  control.effect(() => {
    el.hidden = !source.access.value.can('insert');
    el.disabled = source.edit.value === undefined;
  });
  return control;
}
