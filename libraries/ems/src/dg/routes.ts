/* The platform's own row addresses, resolved in u2 (phase 3 ruling R1): `/domains/<schema>/<table>`
   and `…/<table>/<keyOrId>` open the same `DomainApp` a package app does, so the Dart router needs
   no u2 in the core bundle — it calls the `domainRoutes` func PowerPack exposes over this. An app
   already open over the table takes the address instead of a second one being docked. */
import type * as DG from 'datagrok-api/dg';
import {Control, div} from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';
import {domains} from './facade.js';
import type {DomainTable} from './facade.js';
import {DomainAddress} from './address.js';
import {DomainApp} from './app.js';
import {DomainErrors} from './errors.js';

/** How a table the registry does not know is refused: the platform's registry answers
 * `'validation'` ("Unknown domain table"), the memory backend `'not-found'`. */
const UNKNOWN_TABLE = ['validation', 'not-found'];

/** The view for a `/domains` address, null where the platform keeps the address: anything that is
 * not a table route (`/domains`, `/domains/<schema>` — the schema gallery and the diagram stay
 * Dart) and a table read that failed. A table or schema the registry does not know opens a
 * not-found view with nothing to act on. An app already open over the table takes
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
  let handle: DomainTable;
  try {
    handle = await domains.table(`${schema}.${table}`);
  } catch (e) {
    if (!UNKNOWN_TABLE.includes(DomainErrors.codeOf(e)))
      return null;
    const said = `No table matches "${schema}.${table}". It may have been deleted, renamed, or you may not ` +
      'have access to it.';
    return appView({name: table, content: new Control(div([said], 'u2-async-empty')), path: rest});
  }
  // the platform's own address is the canonical one, and a list opened from it follows the server
  // like the app a package mounts does — two addresses onto one table must not behave differently
  const view = handle.app({path: base, live: true});
  await DomainApp.of(view)!.open(rest);
  return view;
}
