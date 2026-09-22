/* Binding authoring end to end (EMS external bindings, authoring WO-A4): the PowerPack function
   opens the u2 "Create domain schema" dialog over the Northwind connection preset to `orders` —
   on the Design step, the identifier proposed from the connection and the schema; a name › NEXT ›
   VALIDATE › CREATE registers a throwaway schema and opens the u2 app over its `orders`; the
   schema is deleted afterwards. Northwind itself is never written. Skips itself without the
   connection or the CreateDomainSchema privilege. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

const CONNECTION = 'NorthwindBinding:PostgresNorthwind';

category('U2: domain authoring', () => {
  const name = `auth_u2_${Date.now().toString(36)}`;
  let connection: DG.DataConnection | null = null;
  let skip: string | null = null;
  let created = false;
  let flagBefore: boolean | undefined;

  const buttonNamed = (text: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
    .find((b) => b.textContent === text);
  const status = () => document.querySelector('.u2-wizard-status')?.textContent ?? '';
  const reason = () => document.querySelector('.u2-wizard-reason')?.textContent ?? '';

  async function until(check: () => boolean, what: string): Promise<void> {
    for (let i = 0; i < 300 && !check(); i++)
      await delay(100);
    if (!check())
      throw new Error(`timed out waiting for ${what}; status "${status()}", reason "${reason()}"`);
  }

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  before(async () => {
    // grok test's fresh client profile boots with the Beta flag off, and the dialog refuses without it
    flagBefore = grok.shell.settings.enableDomainDatabases;
    grok.shell.settings.enableDomainDatabases = true;
    connection = (await grok.dapi.connections.list({pageSize: 500})).find((c) => c.nqName === CONNECTION) ?? null;
    if (connection === null) {
      skip = `the ${CONNECTION} connection is not on this stand`;
      return;
    }
    try {
      await grok.dapi.domains.draft({connection: CONNECTION, schema: 'public', tables: ['orders']});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'unknown-connection' || e.code === 'forbidden'))
        skip = `${e.code}: ${e.message}`;
      else
        throw e;
    }
  });

  after(async () => {
    for (const view of [...grok.shell.views].filter((v) => (v.path ?? '').includes(`/domains/${name}/`)))
      view.close();
    for (const dialog of document.querySelectorAll<HTMLButtonElement>('.u2-dialog .u2-dialog-close'))
      dialog.click();
    if (created)
      await grok.dapi.domains.schema(name).delete();
    grok.shell.settings.enableDomainDatabases = flagBefore ?? false;
  });

  test('the PowerPack function: connection › design › review creates the schema and opens the app', async () => {
    if (skipped())
      return;
    const running = grok.functions.call('PowerPack:createDomainBinding',
      {connection, schema: 'public', table: 'orders'}) as Promise<string | null>;
    await until(() => /bindable/.test(document.querySelector('.u2-binding-facts')?.textContent ?? ''), 'the draft');
    await until(() => document.querySelector('.u2-manifest-editor') !== null, 'the Design step, opened on');
    expect(document.querySelector('.u2-wizard-step-current')?.textContent, '2Design', 'connection and schema preset');
    const nameInput = document.querySelector<HTMLInputElement>('.u2-manifest-editor [data-u2-name="name"] input')!;
    expect(/^[a-z][a-z0-9_]*$/.test(nameInput.value), true, `a proposed identifier: "${nameInput.value}"`);
    expect(document.querySelector('.u2-manifest-editor [data-u2-name="name"] .u2-input-postfix')?.textContent,
      `registered as ext_${nameInput.value}`);
    nameInput.value = name;
    nameInput.dispatchEvent(new Event('change', {bubbles: true}));
    await until(() => buttonNamed('NEXT')?.disabled === false, 'NEXT after the name');
    buttonNamed('NEXT')!.click();
    await until(() => (document.querySelector('.u2-binding-json')?.textContent ?? '').includes(`"name": "${name}"`),
      'the Review step');
    expect(buttonNamed('CREATE')!.disabled, true, 'CREATE waits for a green validate');
    buttonNamed('VALIDATE')!.click();
    await until(() => status() === 'Validated', 'the dry run');
    expect(buttonNamed('CREATE')!.disabled, false);
    created = true;
    buttonNamed('CREATE')!.click();
    const result = await running;
    expect(result, name);
    expect(document.querySelector('.u2-binding-dialog'), null, 'the dialog closed');
    // not `filter('name = …')`: a schema registered with a friendly name is not found by it
    const schemas = await grok.dapi.domains.schemas.list({pageSize: 500});
    expect(schemas.some((s) => s.name === name), true, 'registered');
    await until(() => [...grok.shell.views].some((v) => (v.path ?? '').includes(`/domains/${name}/orders`)),
      'the u2 app over orders');
    // the fixture's Northwind (the local demo container) ships 830 orders
    expect(await grok.dapi.domains.table(`${name}.orders`).count(), 830, 'Northwind orders through the binding');
  });
});
