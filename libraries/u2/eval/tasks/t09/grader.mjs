import * as grok from 'datagrok-api/grok';
import {clickOk, errors, expect, flush, liveScopes, named, need} from '../../lib/harness.mjs';
import {frame, tableView} from '../../lib/platform.mjs';

async function install(mod, refresh = async () => {}) {
  const table = frame([{name: 'Molecule', type: 'string', semType: 'Molecule'}, {name: 'Stage', type: 'string'}],
    [{Molecule: 'CCO', Stage: 'Idea'}], 'KHL-3');
  const view = tableView(table);
  const session = {refreshes: 0, refresh: () => {
    session.refreshes++;
    return refresh();
  }};
  need(mod, 'campaignRibbon')(view, session);
  await flush();
  const icon = (re, what) => {
    const found = named(view.dart.panels, re);
    expect(found !== null, `no ${what} on the view's ribbon`);
    return found;
  };
  return {table, view, session, icon};
}

const streams = () => grok.events.onViewRemoved.count + grok.events.onCurrentViewChanged.count;

export const checks = {
  'add-row': async (mod) => {
    const {table, view, icon} = await install(mod);
    await clickOk(icon(/add new row/i, 'Add new row'), 'Add new row');
    expect(table.rowCount === 2, `the table has ${table.rowCount} rows after Add new row`);
    expect(errors().length === 0, `Add new row logged: ${errors()[0]}`);
    view.close();
  },
  'refresh-runs': async (mod) => {
    const {session, view, icon} = await install(mod);
    await clickOk(icon(/refresh/i, 'Refresh Campaign'), 'Refresh Campaign');
    expect(session.refreshes === 1, `session.refresh() ran ${session.refreshes} times`);
    view.close();
  },
  'refresh-failure-handled': async (mod) => {
    const {view, icon} = await install(mod, () => Promise.reject(new Error('stand offline')));
    await clickOk(icon(/refresh/i, 'Refresh Campaign'), 'Refresh Campaign');
    await flush();
    const dropped = errors().filter((e) => e.startsWith('unhandled rejection'));
    expect(dropped.length === 0,
      'a failed refresh became an unhandled promise rejection: the action dropped its promise');
    view.close();
  },
  'gone-with-view': async (mod) => {
    const base = liveScopes();
    const subs = streams();
    const {view} = await install(mod);
    view.close();
    await flush();
    expect(liveScopes() === base, `${liveScopes() - base} scopes outlive the closed view`);
    expect(streams() === subs, `${streams() - subs} platform subscriptions outlive the closed view`);
  },
};
