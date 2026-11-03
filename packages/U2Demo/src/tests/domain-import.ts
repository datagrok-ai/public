/* The import wizard against a real domain table (WO 3-9): a frame picked, its columns auto-mapped
   against `apitests.item`, the blocking problems, and Import posting the mapped columns only
   through `DomainTableClient.batch`. Fixture rows are prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {Control, Scope} from '@datagrok-libraries/u2';
import type {Input} from '@datagrok-libraries/u2';
import {domains, DomainTable} from '@datagrok-libraries/ems/src/dg';

const ITEMS = 'apitests.item';

category('U2: domain import', () => {
  const items = () => grok.dapi.domains.table(ITEMS);
  const prefix = `U2I-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  let table: DomainTable;
  /** What the tests build after an await is the category's. */
  let owner: Scope;

  const named = (name: string) =>
    Control.forElement(document.querySelector(`[data-u2-name="${name}"]`)) as Input<string | null> | undefined;
  const buttonNamed = (text: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
    .find((b) => b.textContent === text)!;
  /** Preview's commit button: "Import N …", plain "Import" where the dry run gives no count. */
  const importButton = () => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
    .find((b) => /^Import( |$)/.test(b.textContent ?? ''))!;
  const reason = () => document.querySelector('.u2-wizard-reason')?.textContent ?? '';
  const report = () => document.querySelector('.u2-domain-import-report')?.textContent ?? '';
  /** The preview's result badge per row — what the server's dry run predicts for it. */
  const verdicts = () => [...document.querySelectorAll('.u2-domain-import-rows .u2-data-table-row')]
    .map((row) => row.querySelector('.u2-badge')?.textContent ?? '');

  /** The source frame: every value a string, the way a CSV arrives. `strayColumn` matches no column. */
  function frame(rows: {sku: string, name: string, quantity: string, strayColumn: string}[]): DG.DataFrame {
    return DG.DataFrame.fromColumns(['sku', 'name', 'quantity', 'strayColumn'].map((column) =>
      DG.Column.fromStrings(column, rows.map((r) => String(r[column as keyof typeof rows[0]])))));
  }

  /** The wizard up over the given frame, with the mapping step on screen. The promise is handed
   * back WRAPPED — awaiting an async function that returned it would chain on it. */
  async function open(source: DG.DataFrame): Promise<{running: Promise<unknown>}> {
    const running = owner.build(() => domains.import(table, {source}));
    await waitUntil(() => named('source') !== undefined, 'the wizard is up', 5000);
    buttonNamed('NEXT').click();
    await waitUntil(() => named('sku') !== undefined, 'the mapping step is up');
    return {running};
  }

  before(async () => {
    owner = Scope.root('u2 domain import tests');
    table = await domains.table(ITEMS);
  });

  after(async () => {
    owner.dispose();
    await items().deleteWhere(query);
  });

  test('the mapping is auto-matched, and a column matching nothing is skipped', async () => {
    const {running} = await open(frame([{sku: `${prefix}-1`, name: 'Alpha', quantity: '1', strayColumn: 'x'}]));
    expect(named('sku')!.value, 'sku');
    expect(named('name')!.value, 'name');
    expect(named('quantity')!.value, 'quantity');
    expect(named('strayColumn')!.value, '(skip)');
    named('sku')!.value = '(skip)';
    buttonNamed('NEXT').click();
    await waitUntil(() => reason().includes('Required column "Sku" is not mapped'), reason);
    buttonNamed('CANCEL').click();
    expect(await running, null);
  });

  test('Import posts the mapped columns and the server\'s report comes back', async () => {
    const rows = [1, 2, 3].map((i) => ({sku: `${prefix}-${i}`, name: `Row ${i}`, quantity: String(i), strayColumn: 'x'}));
    const source = frame(rows);
    const {running} = await open(source);
    buttonNamed('NEXT').click();
    await waitUntil(() => verdicts().length === 3 && importButton()?.disabled === false, 'the preview is ready');
    importButton().click();
    await waitUntil(() => report().includes('added'), report);
    expect(report().includes('3 added'), true, report());
    await waitUntil(() => report().includes('Added a "Sku'), report);
    const names = source.columns.names();
    expect(names.slice(0, 4).join(), 'sku,name,quantity,strayColumn', 'the source\'s own columns stay');
    expect(source.col('sku')!.get(0), `${prefix}-1`, 'and so do its values');
    expect(names.length, 5);
    expect(report().includes(`Added a "${names[4]}" column`), true, report());
    expect(source.col(names[4])!.toList().join(), rows.map((r) => r.sku).join(),
      'the column written back holds each row\'s Sku, not its name');
    buttonNamed('CLOSE').click();
    const result = await running as {inserted: number};
    expect(result.inserted, 3);

    const landed = await items().query({filter: query, sort: 'sku'});
    expect(landed.length, 3);
    expect(landed.map((r) => `${r.name}:${r.quantity}`).join(), 'Row 1:1,Row 2:2,Row 3:3',
      'the mapped columns landed, the text coerced by the server');
    expect('strayColumn' in landed[0], false, 'a skipped source column never reached the server');
  });

  test('the preview verdicts are the server\'s dry run, and its counts are the commit\'s', async () => {
    const rows = [1, 2].map((i) =>
      ({sku: `${prefix}-p${i}`, name: `Preview ${i}`, quantity: String(i), strayColumn: 'x'}));
    const {running} = await open(frame(rows));
    buttonNamed('NEXT').click();
    await waitUntil(() => verdicts().length > 0, 'the preview shows its verdicts');
    expect(verdicts().join(), 'Add,Add', 'the result column is what the dry run predicts, as a verb');
    const predicted = verdicts().filter((v) => v === 'Add').length;
    importButton().click();
    await waitUntil(() => report().includes('added'), report);
    expect(report().includes(`${predicted} added`), true, report());
    buttonNamed('CLOSE').click();
    await running;
  });
});
