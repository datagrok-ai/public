/**
 * Resolving the {@link DG.DomainObjectHandler} an app should act through.
 *
 * @module handler
 */

import * as DG from 'datagrok-api/dg';

/**
 * The domain handler that WINS dispatch for [table]: a plugin's
 * {@link DG.DomainObjectHandler} subclass when one is registered for these rows,
 * the reflective platform default otherwise. Every domain member an app calls —
 * actions, detail tabs, editors, openers — goes through it, so an app inherits a
 * plugin's customizations without asking for them, and a table nobody wrote a
 * handler for still works.
 *
 * Rendering is deliberately NOT resolved here: `ui.render` / `ui.renderCard`
 * already dispatch through the platform, which honours handlers of any kind
 * (including non-domain ones) and falls back to the per-table Dart meta.
 */
export function domainHandler(table: string): DG.DomainObjectHandler {
  const own = new DG.DomainObjectHandler(table);
  const resolved = DG.ObjectHandler.forEntity(own.newRow());
  return resolved instanceof DG.DomainObjectHandler && resolved.table === table ? resolved : own;
}

/**
 * Icons of the row actions that leave the row itself alone: Open, History, Copy
 * link and Watch (a subscription of the current user, not row data).
 *
 * The icon — not the name — is the discriminant: a {@link DG.DomainAction}
 * carries no id, and its `name` is a CAPTION that changes with the state
 * ('Watch' ⇄ 'Unwatch'), so keying on it makes two pages disagree about the same
 * action. These are the platform ribbon's own icons.
 */
const READ_ONLY_ACTION_ICONS = ['folder-open', 'history', 'link', 'bell'];

/**
 * Whether running [action] may have changed the row — THE predicate every list
 * and page reloads on, so that Watch and History never pop an unsaved-changes
 * prompt in one place and stay silent in another.
 *
 * An action nobody recognizes counts as changing: reloading needlessly is a
 * round trip, missing a change is stale data.
 */
export function actionChangesRow(action: DG.DomainAction): boolean {
  return !READ_ONLY_ACTION_ICONS.includes(`${action.icon ?? ''}`);
}

/** Whether [action] is the platform's Delete — the one that leaves the page it
 * ran on with nothing to show. Keyed on the icon, see {@link actionChangesRow}. */
export function isDeleteAction(action: DG.DomainAction): boolean {
  return action.icon === 'trash-alt';
}
