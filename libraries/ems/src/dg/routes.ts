/* The platform's own row addresses, resolved in u2 (phase 3 ruling R1): `/domains/<schema>/<table>`
   and `…/<table>/<keyOrId>` open the same `DomainApp` a package app does, so the Dart router needs
   no u2 in the core bundle — it calls the `domainRoutes` func PowerPack exposes over this. An app
   already open over the table takes the address instead of a second one being docked. */
import type * as DG from 'datagrok-api/dg';
import {domains} from './facade.js';
import {DomainAddress} from './address.js';
import {DomainApp} from './app.js';
import {DomainErrors} from './errors.js';

/** The view for a `/domains` address, null where the platform keeps the address: anything that is
 * not a table route (`/domains`, `/domains/<schema>` — the schema gallery and the diagram stay
 * Dart) and a table the registry does not answer for. An app already open over the table takes
 * the address and is answered AS ITSELF — null there would tell the platform no u2 route exists,
 * and the Dart domain view would open over the live app (Back off `?trash=1`). */
export async function route(address: string): Promise<DG.ViewBase | null> {
  const at = address.indexOf('?');
  const match = DomainAddress.ROUTE.exec(at < 0 ? address : address.slice(0, at));
  if (match === null)
    return null;
  const [, schema, table, segment] = match;
  const base = `/domains/${schema}/${table}`;
  const rest = `${base}${segment ?? ''}${at < 0 ? '' : address.slice(at)}`;
  const open = DomainApp.openAt(DomainApp.baseOf(`${schema}.${table}`), rest);
  if (open !== undefined)
    return open;
  // a table the registry refuses, or an app's configuration that throws: the platform's own view
  // still opens, but the reason is not lost — an unknown table is no failure
  const handle = await domains.table(`${schema}.${table}`).catch((e) => {
    if (!DomainErrors.isNotFound(e))
      console.error(`EMS: ${schema}.${table} did not open as a u2 view:`, e);
    return null;
  });
  if (handle === null)
    return null;
  // the platform's own address is the canonical one, and a list opened from it follows the server
  // like the app a package mounts does — two addresses onto one table must not behave differently
  const view = handle.app({path: base, live: true});
  await DomainApp.of(view)!.open(rest);
  return view;
}
