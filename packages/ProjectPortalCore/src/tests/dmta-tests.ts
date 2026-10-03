import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {describeCycle, DmtaEvent, dmtaRows, eventsFromRows, proposalCompound, seriesCycles, testedColumns, testedEvents,
  validateDmta} from '../model/dmta';
import {dmtaOf, registerDmtaSource, unregisterDmtaSource} from '../dmta';
import {loadExampleDesignSets} from '../design-sets';
import {filterValue, ROW_LIMIT} from '../rows';
import {_package} from '../package-test';

const ev = (compound: string, stage: DmtaEvent['stage'], on: string | null, series: string | null = 'S1'): DmtaEvent =>
  ({compound, stage, on, source: 'test', series});

category('Portal: DMTA rules', () => {
  test('one row per compound: each stage reached, its earliest date, the furthest stage first', async () => {
    const rows = dmtaRows([
      ev('C1', 'designed', '2026-04-20T00:00:00Z'), ev('C1', 'requested', '2026-04-22T00:00:00Z'), ev('C1', 'made', null),
      ev('C1', 'tested', '2026-05-18T00:00:00Z'), ev('C1', 'tested', '2026-05-04T00:00:00Z'),
      ev('D7', 'designed', '2026-09-01T00:00:00Z'),
      ev('C2', 'made', null, null), ev('C2', 'analyzed', '2026-06-01T00:00:00Z', null),
    ], (c) => c === 'C2' ? 'S2' : null);
    expect(rows.map((r) => `${r.compound}:${r.stage}`).join(), 'C2:analyzed,C1:tested,D7:designed');
    const c1 = rows.find((r) => r.compound === 'C1')!;
    expect(`${c1.stages.tested} ${c1.stages.made === ''} ${c1.cycleDays}`, '2026-05-04T00:00:00Z true 14');
    expect(rows[0].series, 'S2', 'the series from the program data when no event says it');
    expect(rows[0].stages.designed === undefined, true, 'a later stage implies nothing about an earlier one');
  });

  test('per series: stages reached and the median cycle designed → tested', async () => {
    const rows = dmtaRows([
      ev('A', 'designed', '2026-01-01T00:00:00Z'), ev('A', 'tested', '2026-01-11T00:00:00Z'),
      ev('B', 'designed', '2026-01-01T00:00:00Z'), ev('B', 'tested', '2026-01-21T00:00:00Z'),
      ev('C', 'designed', '2026-01-01T00:00:00Z'), ev('C', 'tested', '2026-01-31T00:00:00Z'),
      ev('D', 'designed', '2026-01-01T00:00:00Z'), ev('E', 'made', null, 'S2')]);
    const [s1, s2] = seriesCycles(rows);
    expect(`${s1.series} ${s1.compounds} ${s1.reached.designed} ${s1.reached.tested} ${s1.medianDays} ${s1.n}`, 'S1 4 4 3 20 3');
    expect(describeCycle(s1), 'S1: 4 compounds · 4 designed, 3 tested · cycle 20 days (median of 3)');
    expect(describeCycle(s2), 'S2: 1 compound · 1 made · no cycle time yet');
  });

  test('sources: a dataset\'s compound and date columns, a proposal\'s compound, a plug-in\'s rows, the config', async () => {
    expect(JSON.stringify(testedColumns(['program', 'compound', 'assay', 'run_date'])), '{"compound":"compound","date":"run_date"}');
    expect(testedColumns(['compound', 'value']), null);
    expect(JSON.stringify(testedColumns(['cpd', 'when'], {compoundColumn: 'cpd', dateColumn: 'when'})), '{"compound":"cpd","date":"when"}');
    const t = testedEvents([{compound: 'C1', run_date: '2026-05-18'}, {compound: 'C1', run_date: '2026-05-04'}, {compound: '', run_date: 'x'}],
      {compound: 'compound', date: 'run_date'}, 'Assays');
    expect(t.map((e) => `${e.compound} ${e.on?.slice(0, 10)}`).join(), 'C1 2026-05-04');
    const regs = [{code: 'EXC-1', smiles: 'CCO'}];
    expect(`${proposalCompound({structure: ' CCO '}, regs)} ${proposalCompound({structure: 'OCC'}, regs)} ` +
      `${proposalCompound({structure: 'C', compound_code: 'EXC-9'}, regs)}`, 'EXC-1 null EXC-9');
    const p = eventsFromRows([{compound: 'C1', stage: 'Requested', on: '2026-02-01'}, {compound: 'C2', stage: 'shipped'}, {stage: 'made'}], 'Req');
    expect(`${p.events.map((e) => `${e.compound}:${e.stage}:${e.source}`).join()} ${p.skipped}`, 'C1:requested:Req 2');
    expect(validateDmta({tested: {dataset: 'x'}, sources: ['nope']}, ['assays']).join('; '),
      'dmta: tested names dataset "x", which the portfolio does not declare; dmta: source "nope" must be a function as Package:Function');
    expect(validateDmta(undefined, []).length, 0);
  });
});

category('Portal: DMTA status', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`, limit: ROW_LIMIT}))[0];
    await loadExampleDesignSets(_package);
  });

  after(async () => {
    unregisterDmtaSource('test requests');
  });

  test('EXS-KIN1: the KIN-A loop from the example design sets, the registry and the assay results', async () => {
    const s = await dmtaOf(smd, kin1);
    const row = (c: string) => s.rows.find((r) => r.compound === c)!;
    expect(`${row('EXC-1001').stage} ${row('EXC-1001').label} ${row('EXC-1001').cycleDays}`, 'tested KIN-A-D01 14');
    expect(row('EXC-1002').label, 'KIN-A-D02', 'matched to the registered compound by its SMILES');
    expect(`${row('KIN-A-D04').stage} ${row('KIN-A-D04').stages.requested?.slice(0, 10)}`, 'requested 2026-04-22',
      'chosen, not registered: requested');
    expect(row('KIN-A-D05').stage, 'designed');
    expect(row('EXC-1101').stage, 'made', 'registered, never designed here, not in the assay results');
    const kinA = s.cycles.find((c) => c.series === 'KIN-A')!;
    expect(`${kinA.medianDays} ${kinA.n}`, '14 3', describeCycle(kinA));
    expect(s.notes.join('; '), '');
  });

  test('a registered source adds requests; the function returns the table', async () => {
    registerDmtaSource({name: 'test requests', events: async (ctx) => ctx.program.code !== 'EXS-KIN1' ? [] :
      [{compound: 'KIN-A-D05', stage: 'requested', on: '2026-04-25T00:00:00Z', source: 'test requests', series: 'KIN-A'}]});
    const s = await dmtaOf(smd, kin1);
    expect(s.rows.find((r) => r.compound === 'KIN-A-D05')!.stage, 'requested');
    unregisterDmtaSource('test requests');
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:dmtaStatus', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
    expect(df.columns.names().slice(0, 9).join(), 'Compound,Design label,Series,Stage,Designed,Requested,Made,Tested,Analyzed');
    expect(df.col('Stage')!.toList().includes('tested'), true);
    const cycles: DG.DataFrame = await grok.functions.call('ProjectPortalCore:dmtaCycles', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
    expect(cycles.col('Series')!.toList().includes('KIN-A'), true);
  });
});
