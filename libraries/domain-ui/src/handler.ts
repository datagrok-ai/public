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
 * A row action with the ONE thing a caller cannot infer: whether running it may
 * have changed the row. Custom actions (the `actions` option of a list or a row
 * page) say so explicitly; the platform's own are classified by icon, see
 * {@link actionChangesRow}.
 */
export interface DomainRowAction extends DG.DomainAction {
  /** Whether running this action may have changed the row — the list and the row
   * page reload on it. Omitted means "decide from the icon". */
  changesRow?: boolean;
}

/** Builds the actions offered for [row]: the handler's permission-gated defaults,
 * passed through [actions] when the caller supplied one (§F.1.6b — append,
 * filter or reorder, without losing the gating). */
export async function rowActions(handler: DG.DomainObjectHandler, row: DG.DomainRow,
  actions?: DomainActionsProvider): Promise<DomainRowAction[]> {
  const defaults = await handler.getRibbonActions(row as any);
  return actions == null ? defaults : actions(row, defaults);
}

/** Customizes the actions offered for a row — see {@link rowActions}. */
export type DomainActionsProvider =
  (row: DG.DomainRow, defaults: DomainRowAction[]) => DomainRowAction[];

/** Icon of the platform's Open action — the one that opens the row's Entity
 * View. */
export const OPEN_ACTION_ICON = 'folder-open';

/**
 * Icons of the row actions that leave the row itself alone: Open, History, Copy
 * link and Watch (a subscription of the current user, not row data).
 *
 * The icon is the FALLBACK discriminant: a platform {@link DG.DomainAction}
 * carries no id, and its `name` is a CAPTION that changes with the state
 * ('Watch' ⇄ 'Unwatch'), so keying on it makes two pages disagree about the same
 * action. These are the platform ribbon's own icons.
 */
const READ_ONLY_ACTION_ICONS = [OPEN_ACTION_ICON, 'history', 'link', 'bell'];

/**
 * Whether running [action] may have changed the row — THE predicate every list
 * and page reloads on, so that Watch and History never pop an unsaved-changes
 * prompt in one place and stay silent in another.
 *
 * An explicit {@link DomainRowAction.changesRow} wins: it is the only way a
 * mutating action with a read-only icon (a custom 'link'-icon action that writes)
 * can be reloaded on. An action nobody recognizes counts as changing: reloading
 * needlessly is a round trip, missing a change is stale data.
 */
export function actionChangesRow(action: DG.DomainAction): boolean {
  const declared = (action as DomainRowAction).changesRow;
  return declared != null ? declared : !READ_ONLY_ACTION_ICONS.includes(`${action.icon ?? ''}`);
}

/** Whether [action] is the platform's Delete — the one that leaves the page it
 * ran on with nothing to show. Keyed on the icon, see {@link actionChangesRow}. */
export function isDeleteAction(action: DG.DomainAction): boolean {
  return action.icon === 'trash-alt';
}

/** Whether [action] opens the row's own page — which the page OF that row has no
 * use for (it is already showing it). */
export function isOpenAction(action: DG.DomainAction): boolean {
  return action.icon === OPEN_ACTION_ICON;
}
