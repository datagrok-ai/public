import {Control, Scope, button, component, computed} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';
import type {DomainSession} from '../core/session.js';
import type {RowView} from '@datagrok-libraries/u2';

/** A session, or a source standing for its own. */
export type SessionTarget = DomainSession | DomainSource;

export interface SessionButtonOptions {
  text?: string;
}

/** A button over a session: a Control rather than a bare button, so the effects it runs have an owner
 * (the view's host disposes ribbon controls with the view). */
export class SessionButton extends Control {
  static readonly WAITING = 'Wait for the batch being saved to finish';

  readonly button: HTMLButtonElement;
  readonly session: DomainSession;

  /** `disabled` says why the button is off, `null` while on — `button()`'s own option, read live. */
  constructor(target: SessionTarget, text: string, options: {primary?: boolean, tooltip: string,
    run: (session: DomainSession) => unknown, disabled: (session: DomainSession) => string | null}) {
    const session = SessionButton.sessionOf(target);
    const scope = new Scope();
    const el = scope.build(() => button(text, () => void options.run(session),
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
    disabled: (session) => session.isSaving.value ? SessionButton.WAITING :
      session.validity.value ?? (session.isDirty.value ? null : 'No changes to save')});
  control.root.dataset.u2 = 'save-button';
  const session = control.session;
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
  const control = new SessionButton(target, options.text ?? 'Discard', {tooltip: 'Discard changes',
    run: (session) => leave?.when() ? leave.run() : session.discard(),
    disabled: (session) => session.isSaving.value ? SessionButton.WAITING :
      session.isDirty.value || (leave?.when() ?? false) ? null : 'Nothing to discard'});
  control.root.dataset.u2 = 'discard-button';
  return control;
}

/** "New": a draft over `values` made current — for a source, the row a create form or a list then
 * edits; a function receives the row that was current (the one just saved, say) so an app carries
 * fields over from one entry to the next. The press is a change of its own — the draft counts
 * against the session from the start, unlike the create page's initial row. Permission ⇒ hidden:
 * absent without the `insert` capability; state ⇒ off, saying "Loading…", until the table is loaded. */
export function newButton(source: DomainSource,
  values: Record<string, unknown> | ((last: RowView | null) => Record<string, unknown>) = {},
  options: SessionButtonOptions = {}): Control {
  const control = component((self) => {
    const el = button(options.text ?? 'New', () => source.newRow(
      typeof values === 'function' ? values(source.currentRow.peek()) : values),
    {disabled: computed(() => source.edit.value === undefined ? 'Loading…' : null)});
    self.effect(() => {
      el.hidden = !source.access.value.can('insert');
    });
    return el;
  });
  control.root.dataset.u2 = 'new-button';
  return control;
}
