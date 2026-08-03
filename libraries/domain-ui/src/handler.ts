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
