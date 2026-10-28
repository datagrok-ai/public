/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {calculateColumns} from '../app/utils/calculate-single-cell';
import {designCampaign, designCsv, expectSame, noComputeTemplate, Sandbox} from './fixtures';

const compute = () => document.querySelector<HTMLElement>('[data-u2-name="compute"]');

/** A compute function of this test alone, registered in the page as `core:<name>`; [run] gets the table it is handed
 * and adds its column under a free name, as Chem's functions do. The molecules come as the column's name: a function
 * registered in the page does not resolve a column from it. */
function testFunction(run: (table: DG.DataFrame) => void): string {
  const name = `zzb3Compute${DG.Utils.randomString(8)}`;
  grok.functions.register({signature: `void ${name}(dataframe table, string molecules)`, run});
  return `core:${name}`;
}

function designTable(): DG.DataFrame {
  const df = DG.DataFrame.fromCsv(designCsv);
  df.col('Molecule')!.semType = DG.SEMTYPE.MOLECULE;
  return df;
}

function computeButton(text: string): HTMLButtonElement | undefined {
  return Array.from(compute()?.querySelectorAll<HTMLButtonElement>('button') ?? [])
    .find((b) => b.textContent?.trim() === text);
}

category('Campaign view: compute', () => {
  test('Calculate additional properties shows every function\'s settings while a table without a name is open', async () => {
    const sb = await Sandbox.create('compute');
    const nameless = DG.DataFrame.fromColumns([DG.Column.fromStrings('a', ['x'])]);
    try {
      const t = noComputeTemplate('ZZB3 Compute', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      grok.shell.addTable(nameless);
      nameless.name = null as unknown as string;
      app.editCompute();
      await waitUntil(() => compute() != null, 'the compute dialog did not open', 10000);
      const rows = Array.from(compute()!.querySelectorAll<HTMLElement>('[role="option"]'));
      expect(rows.length > 1, true, 'the dialog lists functions');
      for (const row of rows)
        row.click();
      expectSame(Array.from(compute()!.querySelectorAll('[role="region"] [data-u2="readonly-field"] .u2-param-source-error'),
        (el) => el.textContent), [], 'the settings that could not be shown');
      computeButton('CANCEL')?.click();
    } finally {
      grok.shell.closeTable(nameless);
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a calculation over the table names its columns after the molecule column and saves no hidden column', async () => {
    const sb = await Sandbox.create('compute');
    try {
      const t = noComputeTemplate('ZZB3 Compute', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      const df = app.dataFrame!;
      // what Chem keeps once a filter or a search read the molecules
      const cache = df.columns.addNewString('~Molecule.canonicalSmiles');
      cache.init((i) => df.get('Molecule', i));
      cache.semType = DG.SEMTYPE.MOLECULE;
      app.editCompute();
      await waitUntil(() => computeButton('OK')?.disabled === false, 'the compute dialog never let OK go', 30000);
      const box = Array.from(compute()!.querySelectorAll('[role="option"]'))
        .find((r) => r.textContent?.trim() === 'Pharmacophores')?.querySelector<HTMLInputElement>('input[type="checkbox"]');
      expect(box != null, true, 'Chem\'s Pharmacophores is offered');
      box!.click();
      computeButton('OK')!.click();
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the calculation was not saved', 60000);
      const calculated = df.columns.names().filter((n) => n.startsWith('Acceptor'));
      expectSame(calculated, ['Acceptor (Molecule)'], 'the calculated column, named after the molecule column');
      expect(app.molColName, 'Molecule', 'the molecule column');
      const saved = await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'));
      expectSame(saved.columns.names().filter((n) => n.includes('~')), [], 'columns of the saved table naming a hidden one');
      const json = await sb.readCampaign('ZZB3-1');
      expectSame([...Object.keys(json.columnSemTypes), ...Object.keys(json.columnTypes), ...Object.keys(json.columnTags)]
        .filter((n) => n.includes('~')), [], 'columns campaign.json names that are hidden ones');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 180000});

  test('a hidden column in a campaign\'s table is dropped on open and is never its molecule column', async () => {
    const sb = await Sandbox.create('compute');
    try {
      const t = noComputeTemplate('ZZB3 Compute', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1,
        columnSemTypes: {'~Molecule.canonicalSmiles': 'Molecule', 'Molecule': 'Molecule'}}),
      '~Molecule.canonicalSmiles,Molecule,Stage,V-iD,Note\nCCO,CCO,Idea,,first\nc1ccccc1,c1ccccc1,Made,,second\n');
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      expect(app.molColName, 'Molecule', 'the molecule column');
      expectSame(app.dataFrame!.columns.names(), ['Molecule', 'Stage', 'V-iD', 'Note'], 'the columns the campaign opens with');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a calculation over the table tells its grid and filters that the values changed', async () => {
    const df = designTable();
    const fn = testFunction((table) => table.columns.addNewInt(table.columns.getUnusedName('ZZB3 length'))
      .init((i) => (table.get('Molecule', i) as string).length));
    await calculateColumns({descriptors: [], externals: {[fn]: {}}}, df, 'Molecule');
    expectSame(df.col('ZZB3 length')?.toList(), [3, 8], 'the first calculation');
    df.set('Molecule', 0, 'CCCO');
    let changed = 0;
    const sub = df.onValuesChanged.subscribe(() => changed++);
    try {
      await calculateColumns({descriptors: [], externals: {[fn]: {}}}, df, 'Molecule');
    } finally {
      sub.unsubscribe();
    }
    expectSame(df.col('ZZB3 length')?.toList(), [4, 8], 'the column filled again in place');
    expect(changed > 0, true, 'the table said its values changed');
  });

  test('the functions read the table\'s other columns, never a hidden one; run again, they fill their own columns', async () => {
    const df = designTable();
    const cache = df.columns.addNewString('~Molecule.canonicalSmiles');
    cache.init((i) => df.get('Molecule', i));
    const seen: string[][] = [];
    const fn = testFunction((table) => {
      seen.push(table.columns.names());
      table.columns.addNewString(table.columns.getUnusedName('Note length')).init((i) => `${table.get('Note', i)}`.length.toString());
    });
    await calculateColumns({descriptors: [], externals: {[fn]: {}}}, df, 'Molecule');
    expectSame(seen[0], ['Molecule', 'Stage', 'V-iD', 'Note'], 'the columns the function is handed');
    expectSame(df.col('Note length')?.toList(), ['5', '6'], 'the column made from the notes');
    df.set('Note', 1, 'the second one');
    await calculateColumns({descriptors: [], externals: {[fn]: {}}}, df, 'Molecule');
    expectSame(df.columns.names(), ['Molecule', 'Stage', 'V-iD', 'Note', '~Molecule.canonicalSmiles', 'Note length'],
      'the columns after the second calculation');
    expectSame(df.col('Note length')?.toList(), ['5', '14'], 'the column filled again from the notes');
  });

  test('fractions calculated for an integer column fill it where it is, its filter and its viewers kept', async () => {
    const df = designTable();
    df.columns.insert(DG.Column.fromInt32Array('MW', new Int32Array([1, 2])), 1);
    const view = grok.shell.addTableView(df);
    try {
      const plot = view.scatterPlot({x: 'MW', y: 'MW'});
      const filters = view.getFiltersGroup({createDefaultFilters: false});
      filters.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'MW', min: 1, max: 2});
      await waitUntil(() => filters.getStates('MW', DG.FILTER_TYPE.HISTOGRAM).length === 1, 'the filter on MW was not added', 5000);
      const fn = testFunction((table) => table.columns.addNewFloat(table.columns.getUnusedName('MW')).init((i) => [46.25, 78.5][i]));
      await calculateColumns({descriptors: [], externals: {[fn]: {}}}, df, 'Molecule');
      expectSame(df.col('MW')?.toList(), [46.25, 78.5], 'the calculated values');
      expect(filters.getStates('MW', DG.FILTER_TYPE.HISTOGRAM).length, 1, 'filters on MW');
      expect(plot.props.xColumnName, 'MW', 'the scatter plot\'s x');
      expectSame(df.columns.names(), ['Molecule', 'MW', 'Stage', 'V-iD', 'Note'], 'the columns, in their order');
    } finally {
      view.close();
      grok.shell.closeTable(df);
    }
  });
});
