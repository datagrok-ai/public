/* Save, Discard and New over a session — a source's own session of one today, the phase-2
   session over several tables tomorrow. State ⇒ disabled: nothing to save, a save in flight, or
   — Save — a batch that may not go, its reason the tooltip; New is the one permission-gated button
   (hidden without `insert`), the fields and actions having degraded by access already. Controls
   rather than bare buttons, so the effects they run have an owner — `appView` disposes ribbon
   controls with the content. */
import {Control, Scope, button, computed} from '@datagrok-libraries/u2';
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

  /** `disabled` is `true`, or why the button is off — `button()`'s own option, read live. */
  constructor(target: SessionTarget, text: string, options: {primary?: boolean, tooltip: string,
    run: (session: DomainSession) => unknown, disabled: (session: DomainSession) => boolean | string}) {
    const session = SessionButton.sessionOf(target);
    const scope = new Scope();
    const el = Scope.runWith(scope, () => button(text, () => void options.run(session),
      {primary: options.primary, tooltip: options.tooltip, disabled: computed(() => options.disabled(session))}));
    super(el, scope);
    this.button = el;
    this.session = session;
  }

  static sessionOf(target: SessionTarget): DomainSession {
    return target instanceof DomainSource ? target.session : target;
  }
}

/** Save: the batch as one transaction through the session, which announces what it saved and
 * hands the focus back to the paired form; the session's primary button, where Tab from a form's
 * last field lands. */
export function saveButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const control = new SessionButton(target, options.text ?? 'Save', {primary: true, tooltip: 'Save (Ctrl+S)',
    run: (session) => session.save(),
    disabled: (session) => session.isSaving.value || (session.validity.value ?? !session.isDirty.value)});
  control.root.dataset.u2 = 'save-button';
  const session = control.session;
  session.primaryButton = control.button;
  control.own(() => {
    if (session.primaryButton === control.button)
      session.primaryButton = undefined;
  });
  return control;
}

export function discardButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const control = new SessionButton(target, options.text ?? 'Discard', {tooltip: 'Discard changes',
    run: (session) => session.discard(), disabled: (session) => !session.isDirty.value || session.isSaving.value});
  control.root.dataset.u2 = 'discard-button';
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
