import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programReadiness, stagesOf, syncPortfolio} from '../portal';
import {COL, portfolioTables, roleColumn, UNSTAGED, unitColumn} from '../portfolio-table';
import {portfolioTableView, VIEWER_TITLES} from '../views/portfolio-table-view';
const portfolioTable = async (code?: string) => (await portfolioTables(code)).programs;
import {_package} from '../package-test';
import {withRestrictedUser} from './restricted-user';

const PORTFOLIOS = ['EX-SMD', 'EX-BIOD', 'EX-BIOC'];
let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  for (const code of PORTFOLIOS)
    await syncPortfolio(code);
})();

const rowOf = (df: DG.DataFrame, program: string, portfolio?: string) => {
  for (let i = 0; i < df.rowCount; i++) {
    if (df.get(COL.program, i) === program && (!portfolio || df.get(COL.portfolio, i) === portfolio))
      return i;
  }
  throw new Error(`${program} is not in ${df.name}`);
};

category('Portal: portfolio table', () => {
  before(async () => { await setUp(); });

  test('one row per program of each example portfolio', async () => {
    for (const code of PORTFOLIOS) {
      const pf = await portfolioByCode(code);
      const programs = await portalDb.programs.query({filter: `portfolio_id = "${pf.id}"`, limit: 10000});
      const df = await portfolioTable(code);
      expect(df.rowCount, programs.length, code);
      expect(df.col(COL.program)!.toList().sort().join(), programs.map((p) => p.code).sort().join(), code);
      expect(df.col(COL.id)!.categories.length, programs.length, `${code}: a program twice`);
      expect(df.col(COL.portfolio)!.categories.join(), code);
    }
  });

  test('the stage is a category in the order of the portfolio\'s model', async () => {
    for (const code of PORTFOLIOS) {
      const pf = await portfolioByCode(code);
      const stages = await stagesOf(pf.stage_model_id);
      const df = await portfolioTable(code);
      const order = [UNSTAGED, ...stages.map((s) => s.code)];
      expect(df.col(COL.stage)!.getCategoryOrder().join(), order.join(), `${code}: stage order`);
      const programs = await portalDb.programs.query({filter: `portfolio_id = "${pf.id}"`, limit: 10000});
      for (const p of programs) {
        const i = rowOf(df, p.code);
        expect(df.get(COL.stage, i), stages.find((s) => s.id === p.stage_id)?.code ?? UNSTAGED, `${code}/${p.code}`);
        expect(df.get(COL.stageOrdinal, i), stages.find((s) => s.id === p.stage_id)?.ordinal ?? 0, `${code}/${p.code}`);
      }
    }
    const dev = await portfolioTable('EX-BIOC');
    expect(dev.get(COL.stage, rowOf(dev, 'EXB-201')), 'D4');
    expect(dev.get(COL.sourceStage, rowOf(dev, 'EXB-201')) !== null, true, 'the registry spelling is kept');
    expect(dev.get(COL.stage, rowOf(dev, 'EXB-205')), UNSTAGED);
  });

  test('readiness agrees with programReadiness', async () => {
    for (const code of PORTFOLIOS) {
      const df = await portfolioTable(code);
      for (let i = 0; i < df.rowCount; i++) {
        const program = df.get(COL.program, i) as string;
        if (df.col(COL.nextStage)!.isNone(i)) {
          expect(df.col(COL.readiness)!.isNone(i), true, `${program} is at its last stage`);
          continue;
        }
        const r = await programReadiness(code, program);
        const at = `${code}/${program} → ${r.stage.code}`;
        expect(df.get(COL.nextStage, i), r.stage.code, at);
        expect(df.get(COL.requiredMet, i), r.readiness.requiredMet, at);
        expect(df.get(COL.requiredTotal, i), r.readiness.requiredTotal, at);
        expect(Math.abs(df.get(COL.readiness, i) - r.readiness.score) < 1e-9, true, at);
        expect(df.get(COL.ready, i), r.readiness.score === 1, at);
      }
    }
  });

  test('the readiness matrix holds the state of every criterion of the next stage', async () => {
    const t = await portfolioTables('EX-SMD');
    const r = await programReadiness('EX-SMD', 'EXS-KIN1');
    const i = rowOf(t.readiness, 'EXS-KIN1');
    for (const c of r.readiness.criteria)
      expect(t.readiness.get(`${r.stage.code}: ${c.code}`, i), c.state, c.code);
    expect(t.readiness.rowCount, t.programs.rowCount);
  });

  test('days in stage count from the last transition into the stage', async () => {
    const df = await portfolioTable('EX-BIOC');
    const pf = await portfolioByCode('EX-BIOC');
    const p = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "EXB-201"`}))[0];
    const log = (await portalDb.stageTransitions.query({filter: `program_id = "${p.id}"`, limit: 10000}))
      .filter((t) => t.to_stage_id === p.stage_id && t.from_stage_id !== t.to_stage_id);
    expect(log.length > 0, true, 'EXB-201 was placed in D4 by a sync');
    const entered = Math.max(...log.map((t) => Number((t.decided_on ?? t.created_on).valueOf())));
    const i = rowOf(df, 'EXB-201');
    expect(df.get(COL.daysInStage, i), Math.floor((Date.now() - entered) / 86400000));
    for (let r = 0; r < df.rowCount; r++)
      expect(df.get(COL.daysInStage, r) >= 0, true, df.get(COL.program, r));
  });

  test('compounds by role, work units and the workspace', async () => {
    const sm = await portfolioTable('EX-SMD');
    const kin = rowOf(sm, 'EXS-KIN1');
    expect(sm.get(roleColumn('lead'), kin), 6);
    const roles = ['explored', 'hit', 'lead', 'backup', 'candidate', 'tool', 'discontinued']
      .reduce((a, r) => a + sm.get(roleColumn(r), kin), 0);
    expect(sm.get(COL.compounds, kin), roles);
    expect(sm.get(unitColumn('campaign'), kin), 2);
    const pf = await portfolioByCode('EX-SMD');
    const p = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "EXS-KIN1"`}))[0];
    expect(sm.get(COL.workspace, kin), !!p.space_id);
    const dev = await portfolioTable('EX-BIOC');
    expect(dev.get(unitColumn('process_stage'), rowOf(dev, 'EXB-201')), 4);
    expect(dev.col(unitColumn('campaign')), null, 'EX-BIOC has no campaigns');
  });

  test('a filter reduces the rows', async () => {
    const df = await portfolioTable('EX-BIOC');
    const d4 = df.col(COL.stage)!.toList().filter((s) => s === 'D4').length;
    expect(d4 > 0 && d4 < df.rowCount, true, `D4: ${d4} of ${df.rowCount}`);
    df.filter.init((i) => df.get(COL.stage, i) === 'D4');
    expect(df.filter.trueCount, d4);
    df.filter.and(DG.BitSet.create(df.rowCount, (i) => df.get(COL.ready, i) === true));
    expect(df.filter.trueCount <= d4, true);
  });

  test('every portfolio in one table: the portfolio column filters it', async () => {
    const all = await portfolioTable();
    for (const code of PORTFOLIOS) {
      const one = await portfolioTable(code);
      all.filter.init((i) => all.get(COL.portfolio, i) === code);
      expect(all.filter.trueCount, one.rowCount, code);
    }
    expect(all.col(COL.stage)!.getCategoryOrder()[0], UNSTAGED);
  });

  test('the table view opens with the filters and its viewers, and they follow the filter', async () => {
    const tv = await portfolioTableView('EX-BIOC');
    try {
      const types = Array.from(tv.viewers).map((v) => v.type);
      for (const t of [DG.VIEWER.GRID, DG.VIEWER.FILTERS, DG.VIEWER.BAR_CHART, DG.VIEWER.BOX_PLOT])
        expect(types.includes(t), true, `${t} in ${types.join()}`);
      const t = (tv as any).portfolioTables;
      const extra = (tv as any).extraViewers as DG.Viewer[];
      expect(extra.map((v) => v.dataFrame.name).join(), [t.readiness.name, t.compounds.name].join());
      expect(tv.path.endsWith('/Portal/EX-BIOC/table'), true, tv.path);
      const df = tv.dataFrame;
      df.filter.init((i) => df.get(COL.stage, i) === 'D4');
      df.filter.fireChanged();
      await DG.delay(200);
      expect(t.readiness.filter.trueCount, df.filter.trueCount, 'the readiness matrix follows the filter');
    } finally {
      tv.close();
    }
  });

  test('a restricted user\'s table holds only the portfolios they are in', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const ran = await withRestrictedUser('pptable', async (u) => {
      expect((await u.asUser(() => portfolioTable())).rowCount, 0, 'an outsider sees programs');
      const viewers = await grok.dapi.groups.find(sm.viewers!);
      const personal = await grok.dapi.groups.find(u.group);
      await grok.dapi.groups.includeTo(personal, viewers);
      try {
        const df = await u.asUser(() => portfolioTable());
        const seen = df.col(COL.portfolio)!.categories.filter((c) => c);
        expect(seen.every((c) => c === 'EX-SMD'), true, `saw ${seen.join()}`);
      } finally {
        await grok.dapi.groups.excludeFrom(personal, viewers);
      }
    });
    expect(ran !== null, true, 'no restricted session could be created on this stand');
  });
});
