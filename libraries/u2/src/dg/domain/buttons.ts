/* Save, Discard and New over a session — a source's own session of one today, the phase-2
   session over several tables tomorrow. State ⇒ disabled: nothing to save, or a save in flight;
   New is the one permission-gated button (hidden without `insert`), the fields and actions
   having degraded by access already. Controls rather than bare buttons, so the effects they run
   have an owner — `appView` disposes ribbon controls with the content. */
import {Control} from '../../core/component.js';
import {button} from '../../core/elements.js';
import {DomainSource} from '../../sources/domain-source.js';
import type {DomainSession} from '../../sources/session.js';
import type {RowView} from '../../sources/rows-like.js';

/** A session, or a source standing for its own. */
export type SessionTarget = DomainSession | DomainSource;

export interface SessionButtonOptions {
  text?: string;
}

export class SessionButton extends Control {
  readonly button: HTMLButtonElement;
  readonly session: DomainSession;

  constructor(target: SessionTarget, text: string,
    options: {primary?: boolean, run: (session: DomainSession) => unknown}) {
    const session = SessionButton.sessionOf(target);
    const el = button(text, () => void options.run(session), {primary: options.primary});
    super(el);
    this.button = el;
    this.session = session;
    this.effect(() => el.disabled = !session.isDirty.value || session.isSaving.value);
  }

  static sessionOf(target: SessionTarget): DomainSession {
    return target instanceof DomainSource ? target.session : target;
  }
}

/** Save: the batch as one transaction through the session, which announces what it saved and
 * hands the focus back to the paired form. */
export function saveButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const control = new SessionButton(target, options.text ?? 'Save', {primary: true, run: (session) => session.save()});
  control.root.dataset.u2 = 'save-button';
  return control;
}

export function discardButton(target: SessionTarget, options: SessionButtonOptions = {}): SessionButton {
  const control = new SessionButton(target, options.text ?? 'Discard', {run: (session) => session.discard()});
  control.root.dataset.u2 = 'discard-button';
  return control;
}

/** "New": a pristine draft over `values` made current — for a source, the row a create form or a
 * list then edits; a function receives the row that was current (the one just saved, say) so an
 * app carries fields over from one entry to the next. Permission ⇒ hidden: absent without the
 * `insert` capability; state ⇒ disabled until the table is loaded. */
export function newButton(source: DomainSource,
  values: Record<string, unknown> | ((last: RowView | null) => Record<string, unknown>) = {},
  options: SessionButtonOptions = {}): Control {
  const el = button(options.text ?? 'New', () => source.newRow(
    typeof values === 'function' ? values(source.currentRow.peek()) : values, {pristine: true}));
  const control = new Control(el);
  control.root.dataset.u2 = 'new-button';
  control.effect(() => {
    el.hidden = !source.access.value.can('insert');
    el.disabled = source.edit.value === undefined;
  });
  return control;
}
