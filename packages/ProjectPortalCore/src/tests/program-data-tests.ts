import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {validatePortfolio} from '../model/config';
import {
  ColumnInfo, datasetTableName, defaultViewers, MAIN_DATASET, programDatasets, validateProgramDatasets,
} from '../model/program-data';
import {datasetsOf, openProgramTable, programTable, viewersFor} from '../program-data';
import {_package} from '../package-test';

const extra = (over: Record<string, unknown> = {}) =>
  ({code: 'assays', name: 'Assays', connection: 'C', sql: 'SELECT 1 WHERE p = @program', ...over}) as any;

category('Portal: program data rules', () => {
  test('EX-SMD names its program dataset and declares assay results beside it', async () => {
    const smd = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-SMD')!;
    expect(validatePortfolio(smd).join('; '), '');
    const all = programDatasets(smd.programData, smd.programDatasets);
    expect(all.map((d) => `${d.code}:${d.name}`).join(), 'program:SAR table,assays:Assay results');
  });

  test('the program dataset comes first and is called Program data unless named; none without SQL', async () => {
    const all = programDatasets({connection: 'C', sql: 'S @program'}, [extra()]);
    expect(all.map((d) => d.name).join(), 'Program data,Assays');
    expect(programDatasets(null, [extra()]).map((d) => d.code).join(), 'assays');
    expect(programDatasets(undefined, undefined).length, 0);
  });

  test('a further dataset needs a unique code, a name, a connection, SQL with @program, known viewers', async () => {
    expect(validateProgramDatasets([extra()]).length, 0);
    const errors = validateProgramDatasets([extra({code: MAIN_DATASET}), extra(), extra({sql: 'SELECT 1'}),
      extra({code: 'x', name: '', viewers: [{type: 'Radar'}]})], [{type: 'Sunburst'} as any]);
    for (const e of ['code is used twice ("program" is the program dataset)', 'programDatasets assays: code is used twice',
      'with @program', 'x: name is required', 'viewer "Radar"', 'programData: viewer "Sunburst"'])
      expect(errors.some((x) => x.includes(e)), true, `${e} in: ${errors.join('; ')}`);
  });

  test('the table is named after the program, and the dataset when it is a further one', async () => {
    expect(datasetTableName({code: 'P1', name: 'Kinase'}, {code: MAIN_DATASET, name: 'SAR table'}), 'P1 · Kinase');
    expect(datasetTableName({code: 'P1', name: 'P1'}, {code: 'assays', name: 'Assays'}), 'P1 · Assays');
  });

  test('default viewers: discovery plots two measures, development counts a category; never the program column', async () => {
    const cols: ColumnInfo[] = [{name: 'program', kind: 'string', categories: 1}, {name: 'compound', kind: 'string', categories: 40},
      {name: 'series', kind: 'string', categories: 3}, {name: 'ic50', kind: 'number'}, {name: 'sol', kind: 'number'}];
    const disc = defaultViewers('discovery', 'small_molecule', cols);
    expect(disc.map((v) => v.type).join(), 'Scatter plot,Histogram');
    expect(JSON.stringify(disc[0].options), '{"xColumnName":"ic50","yColumnName":"sol","colorColumnName":"series"}');
    const dev = defaultViewers('development', 'biologic', [{name: 'program', kind: 'string', categories: 2},
      {name: 'purpose', kind: 'string', categories: 3}, {name: 'released', kind: 'other'}]);
    expect(dev.map((v) => `${v.type}:${JSON.stringify(v.options)}`).join(), 'Bar chart:{"splitColumnName":"purpose"}');
    expect(defaultViewers('discovery', null, [{name: 'name', kind: 'string', categories: 5}]).length, 0);
  });
});

let kin: ProgramRow;
const opened: DG.ViewBase[] = [];

category('Portal: program data', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    const smd = await portfolioByCode('EX-SMD');
    kin = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
    expect(!!kin, true, 'EXS-KIN1 is synced (syncPortfolio EX-SMD)');
  });

  after(async () => {
    for (const v of opened.splice(0))
      v.close();
  });

  test('the program dataset of EXS-KIN1: its rows only, named after it, tagged', async () => {
    const t = await programTable('EX-SMD', 'EXS-KIN1');
    expect(t.df.rowCount, 11);
    expect(t.df.name, datasetTableName(kin, {code: MAIN_DATASET, name: 'SAR table'}));
    expect(t.df.col('program')!.categories.join(), 'EXS-KIN1');
    expect(`${t.df.getTag('portal.portfolio')} ${t.df.getTag('portal.program')} ${t.df.getTag('portal.dataset')}`,
      'EX-SMD EXS-KIN1 program');
    expect(viewersFor(t).map((v) => v.type).join(), 'Scatter plot,Histogram', 'discovery defaults');
  });

  test('a further dataset by code; an unknown one and a portfolio without data are refused', async () => {
    const t = await programTable('EX-SMD', 'EXS-KIN1', 'assays');
    expect(t.df.rowCount, 8);
    expect(t.df.name.endsWith(' · Assay results'), true, t.df.name);
    expect(viewersFor(t).map((v) => v.type).join(), 'Box plot', 'the viewers it declares');
    const refusal = async (f: () => Promise<unknown>) => {
      try {
        await f();
        return '';
      } catch (e: any) {
        return String(e?.message ?? e);
      }
    };
    expect(await refusal(() => programTable('EX-SMD', 'EXS-KIN1', 'nope')), 'EX-SMD has no dataset nope; it has program, assays');
    expect(await refusal(() => programTable('EX-BIOD', 'ANY')), 'EX-BIOD has no program ANY, or it is not shared with you');
    const biod = await portfolioByCode('EX-BIOD');
    expect(datasetsOf(biod).length, 0, 'EX-BIOD declares no program data');
  });

  test('Open data opens a table view named after the program, with its viewers', async () => {
    const tv = await openProgramTable('EX-SMD', 'EXS-KIN1');
    opened.push(tv);
    expect(grok.shell.tv?.dataFrame?.name, datasetTableName(kin, {code: MAIN_DATASET, name: 'SAR table'}));
    const types = Array.from(tv.viewers).map((v) => v.type);
    expect(types.includes(DG.VIEWER.SCATTER_PLOT) && types.includes(DG.VIEWER.HISTOGRAM), true, types.join());
  });
});
