/* The import wizard against a real domain table (WO 3-9): a frame picked, its columns auto-mapped
   against `apitests.item`, the blocking problems, and FINISH posting the mapped columns only
   through `DomainTableClient.batch`. Fixture rows are prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Control} from '@datagrok-libraries/u2';
import {domains, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const ITEMS = 'apitests.item';

category('U2: domain import', () => {
  const items = () => grok.dapi.domains.table(ITEMS);
  const prefix = `U2I-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  let table: DomainTable;

  const named = (name: string) =>
    Control.forElement(document.querySelector(`[data-u2-name="${name}"]`)) as unknown as
      {value: {value: string | null}} | undefined;
  const buttonNamed = (text: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
    .find((b) => b.textContent === text)!;
  const reason = () => document.querySelector('.u2-wizard-reason')?.textContent ?? '';
  const report = () => document.querySelector('.u2-domain-import-report')?.textContent ?? '';
  /** The preview's leading cell per row — what the server's dry run predicts for it. */
  const verdicts = () => [...document.querySelectorAll('.u2-domain-import-rows .u2-data-table-row')]
    .map((row) => row.firstElementChild?.textContent ?? '');

  /** The source frame: every value a string, the way a CSV arrives. `strayColumn` matches no column. */
  function frame(rows: {sku: string, name: string, quantity: string, strayColumn: string}[]): DG.DataFrame {
    return DG.DataFrame.fromColumns(['sku', 'name', 'quantity', 'strayColumn'].map((column) =>
      DG.Column.fromStrings(column, rows.map((r) => String(r[column as keyof typeof rows[0]])))));
  }

  /** The wizard up over the given frame, with the mapping step on screen. The promise is handed
   * back WRAPPED — awaiting an async function that returned it would chain on it. */
  async function open(source: DG.DataFrame): Promise<{running: Promise<unknown>}> {
    const running = domains.import(table, {source});
    for (let i = 0; i < 100 && named('source') === undefined; i++)
      await delay(50);
    expect(named('source') !== undefined, true, 'the wizard is up');
    buttonNamed('NEXT').click();
    await delay(200);
    return {running};
  }

  before(async () => {
    table = await domains.table(ITEMS);
  });

  after(async () => {
    await items().deleteWhere(query);
  });

  test('the mapping is auto-matched, and a column matching nothing is skipped', async () => {
    const {running} = await open(frame([{sku: `${prefix}-1`, name: 'Alpha', quantity: '1', strayColumn: 'x'}]));
    expect(named('sku')!.value.value, 'sku');
    expect(named('name')!.value.value, 'name');
    expect(named('quantity')!.value.value, 'quantity');
    expect(named('strayColumn')!.value.value, '(skip)');
    named('sku')!.value.value = '(skip)';
    await delay(100);
    buttonNamed('NEXT').click();
    await delay(100);
    expect(reason().includes('Required column "Sku" is not mapped'), true, reason());
    buttonNamed('CANCEL').click();
    expect(await running, null);
  });

  test('FINISH posts the mapped columns and the server\'s report comes back', async () => {
    const rows = [1, 2, 3].map((i) => ({sku: `${prefix}-${i}`, name: `Row ${i}`, quantity: String(i), strayColumn: 'x'}));
    const {running} = await open(frame(rows));
    buttonNamed('NEXT').click();
    await delay(400);
    buttonNamed('NEXT').click();
    for (let i = 0; i < 100 && !report().includes('inserted'); i++)
      await delay(50);
    expect(report().includes('3 inserted, 0 updated'), true, report());
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
    for (let i = 0; i < 100 && verdicts().length === 0; i++)
      await delay(50);
    expect(verdicts().join(), 'Add,Add', 'the leading column is what the dry run predicts, as a verb');
    const predicted = verdicts().filter((v) => v === 'Add').length;
    buttonNamed('NEXT').click();
    for (let i = 0; i < 100 && !report().includes('inserted'); i++)
      await delay(50);
    expect(report().includes(`${predicted} inserted, 0 updated`), true, report());
    buttonNamed('CLOSE').click();
    await running;
  });
});
