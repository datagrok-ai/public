/* The platform's own row addresses, resolved in u2 (phase 3 ruling R1): `/domains/<schema>/<table>`
   and `…/<table>/<keyOrId>` open the same `DomainApp` a package app does, so the Dart router needs
   no u2 in the core bundle — it calls the `domainRoutes` func PowerPack exposes over this. An app
   already open over the table takes the address instead of a second one being docked. */
import type * as DG from 'datagrok-api/dg';
import {domains} from './index.js';
import {DomainApp} from './app.js';

/** `/domains/<schema>/<table>[/<keyOrId>]`; the query is `DomainApp.open`'s to read. */
const ADDRESS = /^\/domains\/(\w+)\/(\w+)(\/[^/?#]*)?$/i;

/** The view for a `/domains` address, null where the platform keeps the address: anything that is
 * not a table route (`/domains`, `/domains/<schema>` — the schema gallery and the diagram stay
 * Dart), a table the registry does not answer for, and an address an app already open took. */
export async function route(address: string): Promise<DG.ViewBase | null> {
  const at = address.indexOf('?');
  const match = ADDRESS.exec(at < 0 ? address : address.slice(0, at));
  if (match === null)
    return null;
  const [, schema, table, segment] = match;
  const base = `/domains/${schema}/${table}`;
  const rest = `${base}${segment ?? ''}${at < 0 ? '' : address.slice(at)}`;
  if (DomainApp.openAt(DomainApp.baseOf(`${schema}.${table}`), rest))
    return null;
  const handle = await domains.table(`${schema}.${table}`).catch(() => null);
  if (handle === null)
    return null;
  const view = handle.app({path: base});
  await DomainApp.of(view)!.open(rest);
  return view;
}
