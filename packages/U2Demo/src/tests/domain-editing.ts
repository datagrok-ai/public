/* Binding editing end to end (EMS live-binding editing): the PowerPack function opens the u2
   "Edit binding" dialog over a throwaway binding registered from the Northwind draft (`orders`
   alone) — on the Design step over the registered manifest, a no-op gated; a candidate table
   included, a column restricted and a group granted are listed on Review, VALIDATE annotates
   them, SAVE lands them in one apply and the snapshot shows it; a binding moved meanwhile makes
   SAVE reload it with the report. Northwind itself is never written: every change is a registry
   row. Skips itself without the connection or the CreateDomainSchema privilege. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Control} from '@datagrok-libraries/u2';
import type {ManifestEditor} from '@datagrok-libraries/u2/src/dg/domain/authoring/manifest-editor.js';

const CONNECTION = 'NorthwindBinding:PostgresNorthwind';

category('U2: domain editing', () => {
  const name = `edit_u2_${Date.now().toString(36)}`;
  const handle = grok.dapi.domains.schema(name);
  let skip: string | null = null;
  let created = false;
  let flagBefore: boolean | undefined;

  // every read is scoped to this wizard's root (`.u2-binding-dialog`, inside the platform dialog):
  // the harness page may hold other u2 dialogs and wizards
  const WIZARD = '.u2-binding-dialog';
  const wizard = () => document.querySelector<HTMLElement>(WIZARD);
  const buttonNamed = (text: string) => [...wizard()?.querySelectorAll<HTMLButtonElement>('button') ?? []]
    .find((b) => b.textContent === text);
  const status = () => wizard()?.querySelector('.u2-wizard-status')?.textContent ?? '';
  const reason = () => wizard()?.querySelector('.u2-wizard-reason')?.textContent ?? '';
  const rows = (within: string) => [...wizard()?.querySelectorAll(`${within} .u2-binding-change`) ?? []]
    .map((r) => [...r.children].map((c) => c.textContent));
  const editor = (): ManifestEditor => Control.forElement(wizard()?.querySelector('.u2-manifest-editor') ?? null) as ManifestEditor;
  const state = () => `wizard: ${wizard() !== null}, step "${wizard()?.querySelector('.u2-wizard-step-current')?.textContent ?? ''}", ` +
    `status "${status()}", reason "${reason()}"`;

  // by the clock, not by polls: a background page's timers are throttled, and the wait must
  // name its state before the harness's own timeout does
  async function until(check: () => boolean, what: string): Promise<void> {
    const started = Date.now();
    let reported = started;
    while (!check()) {
      const now = Date.now();
      if (now - started > 25000)
        throw new Error(`timed out waiting for ${what} — ${state()}`);
      if (now - reported > 3000) {
        reported = now;
        console.log(`domain editing: still waiting for ${what} (${((now - started) / 1000).toFixed(0)} s) — ${state()}`);
      }
      await delay(100);
    }
    console.log(`domain editing: ${what} (${((Date.now() - started) / 1000).toFixed(1)} s)`);
  }

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  const tokens = async (): Promise<{ifVersion: string; ifIncarnation: string}> => {
    const m = await handle.manifest();
    return {ifVersion: m.version, ifIncarnation: m.incarnation};
  };

  /** Opens the dialog and waits for Design; `running` is the function's own promise, handed back
   * inside an object — returned bare, an async function would await it, i.e. the close. */
  async function open(): Promise<{running: Promise<boolean>}> {
    const open = [...document.querySelectorAll<HTMLElement>('.u2-dialog')];
    console.log(`domain editing: opening the dialog over ${name}; u2 dialogs open before: ${open.length} ` +
      `[${open.map((d) => d.querySelector('.u2-dialog-title')?.textContent ?? '?').join(' | ')}]`);
    for (const close of document.querySelectorAll<HTMLButtonElement>('.u2-dialog .u2-dialog-close'))
      close.click();
    const running = grok.functions.call('PowerPack:editDomainBinding', {schema: name}) as Promise<boolean>;
    await until(() => wizard()?.querySelector('.u2-manifest-editor') != null, 'the Design step');
    const e = editor();
    console.log(`domain editing: catalog ${e.model.catalog}, ${e.model.blockers.value.length} blockers; ${state()}`);
    expect(wizard()!.querySelector('.u2-wizard-step-current')?.textContent, '2Design', 'opened on Design');
    return {running};
  }

  before(async () => {
    // grok test's fresh client profile boots with the Beta flag off, and the dialog refuses without it
    flagBefore = grok.shell.settings.enableDomainDatabases;
    grok.shell.settings.enableDomainDatabases = true;
    if (!(await grok.dapi.connections.list({pageSize: 500})).some((c) => c.nqName === CONNECTION)) {
      skip = `the ${CONNECTION} connection is not on this stand`;
      return;
    }
    let draft: DG.DomainDraft;
    try {
      draft = await grok.dapi.domains.draft({connection: CONNECTION, schema: 'public', tables: ['orders']});
      console.log('domain editing: the Northwind draft is in');
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Editing probe', manifest: draft.manifest});
      created = true;
      console.log(`domain editing: ${name} registered`);
    } catch (e: any) {
      if (e instanceof DG.DomainError && ['unknown-connection', 'forbidden'].includes(e.code))
        skip = `${e.code}: ${e.message}`;
      else
        throw e;
    }
  });

  after(async () => {
    for (const view of [...grok.shell.views].filter((v) => (v.path ?? '').includes(`/domains/${name}/`)))
      view.close();
    for (const close of document.querySelectorAll<HTMLButtonElement>('.u2-dialog .u2-dialog-close'))
      close.click();
    if (created)
      await handle.delete();
    grok.shell.settings.enableDomainDatabases = flagBefore ?? false;
  });

  test('the PowerPack function: a no-op is gated; a candidate table, a restriction and a grant reviewed, validated and saved in one apply', async () => {
    if (skipped())
      return;
    const {running} = await open();
    await until(() => reason() === 'Nothing changed', 'the no-op gate');
    expect(buttonNamed('NEXT')!.disabled, true);
    const e = editor();
    expect(e.model.editing, true);
    expect(e.model.friendlyName.value, 'Editing probe', 'the entity\'s caption');
    expect(e.model.tables.value.filter((t) => t.included).map((t) => t.remote).join(','), 'orders');
    const allUsers = await grok.dapi.groups.filter('friendlyName = "All users"').first();
    e.model.includeTable('shippers', true);
    e.access.addGroup({kind: 'table', table: 'orders'}, {id: allUsers.id, label: allUsers.friendlyName});
    e.access.setVisibility('orders', 'freight', []);
    await until(() => buttonNamed('NEXT')?.disabled === false, 'NEXT after the edits');
    buttonNamed('NEXT')!.click();
    await until(() => wizard()?.querySelector('.u2-binding-changes') != null, 'the Review step');
    const listed = rows('.u2-binding-changes').map((r) => r[0]);
    expect(listed.slice(0, 2).join('|'), 'Table shippers added|orders: View for All users');
    // the author is kept on the column they restrict, with Edit too where they may edit the table
    expect(/^orders\.freight restricted — visible to you alone/.test(listed[2] ?? ''), true, listed.join('|'));
    expect(buttonNamed('SAVE')!.disabled, true, 'SAVE waits for a green validate');
    buttonNamed('VALIDATE')!.click();
    await until(() => status() === 'Validated', 'the dry run');
    expect(buttonNamed('SAVE')!.disabled, false);
    const before = await tokens();
    buttonNamed('SAVE')!.click();
    await until(() => wizard()?.querySelector('.u2-binding-saved') != null, 'the Saved step');
    expect(wizard()!.querySelector('.u2-binding-created-head')?.textContent,
      `SavedDomain schema ${name} is at version ${Number(before.ifVersion) + 1}`);
    const [manifest, access] = await Promise.all([handle.manifest(), handle.access()]);
    expect(Object.keys(manifest.tables).sort().join(','), 'orders,shippers', 'the candidate is registered');
    expect(access.tables.orders.grants!.some((g) => g.group.id === allUsers.id && g.permissions.includes('View')), true,
      JSON.stringify(access.tables.orders.grants));
    expect(access.columns['orders.freight'].state, 'restricted');
    buttonNamed('CLOSE')!.click();
    expect(await running, true);
    expect(wizard(), null, 'the dialog closed');
  }, {timeout: 120000});

  test('a binding moved meanwhile: SAVE reloads it and reports what was kept, what conflicts, what is gone', async () => {
    if (skipped())
      return;
    const {running} = await open();
    await until(() => reason() === 'Nothing changed', 'the no-op gate');
    const e = editor();
    e.model.setDescription('edited in the dialog');
    e.model.setFriendlyName('orders', 'Orders (dialog)');
    await until(() => buttonNamed('NEXT')?.disabled === false, 'NEXT after the edits');
    buttonNamed('NEXT')!.click();
    await until(() => wizard()?.querySelector('.u2-binding-changes') != null, 'the Review step');
    buttonNamed('VALIDATE')!.click();
    await until(() => status() === 'Validated', 'the dry run');
    // the same caption changed behind the dialog: a conflict; the description is the dialog's alone
    const registered = await handle.manifest();
    const t = {ifVersion: registered.version, ifIncarnation: registered.incarnation};
    await handle.apply({...t, tables: {orders: {...registered.tables.orders, friendlyName: 'Orders (elsewhere)'}}});
    buttonNamed('SAVE')!.click();
    await until(() => /^Reloaded at version/.test(status()), 'the reload');
    expect(status(), `Reloaded at version ${Number(t.ifVersion) + 1}: 1 edit kept, 1 conflict, 0 dropped — validate again`);
    expect(rows('.u2-binding-reloaded').map((r) => r.join(' — ')).join('|'),
      'Table orders: friendly name "Orders (dialog)" — the server\'s "Orders (elsewhere)" stands; yours was "Orders (dialog)"');
    expect(reason(), 'Look over the reloaded edits on Design');
    expect(e.model.table('orders')!.friendlyName, 'Orders (elsewhere)', 'the server\'s value stands');
    buttonNamed('CANCEL')!.click();
    expect(await running, false);
  }, {timeout: 120000});
});
