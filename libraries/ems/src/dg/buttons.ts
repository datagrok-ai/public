/* Save, Discard and New over a session — a source's own session of one today, the phase-2
   session over several tables tomorrow. State ⇒ disabled: nothing to save, or a save in flight;
   Save also while the session's validity blocks it (a field in error, a refusal not yet answered);
   New is the one permission-gated button (hidden without `insert`), the fields and actions
   having degraded by access already. Controls rather than bare buttons, so the effects they run
   have an owner — `appView` disposes ribbon controls with the content. */
import {Control, button, iconButton} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';
import type {DomainSession} from '../core/session.js';
import type {RowView} from '@datagrok-libraries/u2';

/** A session, or a source standing for its own. */
export type SessionTarget = DomainSession | DomainSource;

export interface SessionButtonOptions {
  text?: string;
  /** An icon button with the text as its tooltip — the platform ribbon's form; text by default. */
  icon?: string;
}

export class SessionButton extends Control {
  readonly button: HTMLButtonElement;
  readonly session: DomainSession;

  constructor(target: SessionTarget, text: string,
    options: {primary?: boolean, icon?: string, gated?: boolean, run: (session: DomainSession) => unknown}) {
    const session = SessionButton.sessionOf(target);
    const run = () => void options.run(session);
    const el = options.icon === undefined ? button(text, run, {primary: options.primary}) :
      iconButton(options.icon, run, {tooltip: text});
    super(el);
    this.button = el;
    this.session = session;
    this.effect(() => el.disabled = !session.isDirty.value || session.isSaving.value ||
      options.gated === true && session.validity.value !== null);
  }

  static sessionOf(target: SessionTarget): DomainSession {
    return target instanceof DomainSource ? target.session : target;
  }
}

/** Save: the batch as one transaction through the session, which announces what it saved and
 * hands the focus back to the paired form; the session's primary button, where Tab from a form's
 * last field lands. */
export function saveButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const control = new SessionButton(target, options.text ?? 'Save',
    {primary: true, gated: true, run: (session) => session.save()});
  control.root.dataset.u2 = 'save-button';
  control.button.title = 'Save (Ctrl+S)';
  const session = control.session;
  session.primaryButton = control.button;
  control.own(() => {
    if (session.primaryButton === control.button)
      session.primaryButton = undefined;
  });
  return control;
}

export function discardButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const icon = options.icon;
  const control = new SessionButton(target, options.text ?? (icon === undefined ? 'Discard' : 'Discard changes'),
    {icon, run: (session) => session.discard()});
  control.root.dataset.u2 = 'discard-button';
  if (icon === undefined)
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
  const add = () => source.newRow(typeof values === 'function' ? values(source.currentRow.peek()) : values);
  const el = options.icon === undefined ? button(options.text ?? 'New', add) :
    iconButton(options.icon, add, {tooltip: options.text ?? 'New'});
  const control = new Control(el);
  control.root.dataset.u2 = 'new-button';
  control.effect(() => {
    el.hidden = !source.access.value.can('insert');
    el.disabled = source.edit.value === undefined;
  });
  return control;
}
