import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {planWorkUnits, validateWorkUnits, WorkUnitPlanInput, WorkUnitRef} from '../model/work-units';
import {portalDb} from '../generated/db';
import {enabledWorkUnits, loadAllPortfolios, portfolioByCode, syncPortfolio, syncWorkUnits, workUnitsOf} from '../portal';
import {_package} from '../package-test';

const programs = new Map([['P1', 'id-1'], ['P2', 'id-2']]);
const plan = (over: Partial<WorkUnitPlanInput>) => planWorkUnits({portfolio: 'X', kind: 'series', feed: 'f',
  rows: [], programs, existing: [], ...over});
const unit = (programId: string, code: string, values: Record<string, unknown>, extra: Partial<WorkUnitRef> = {}):
  WorkUnitRef => ({id: `u-${code}`, programId, code, source: 'feed', missingFromFeed: false,
  values: {name: code, details: null, ...values}, ...extra});

category('Portal: work-unit rules', () => {
  test('a portfolio enables the units it uses, and feeds only those', async () => {
    expect(validateWorkUnits(['campaign', 'series'], [], ['hit', 'lead']).join(), '');
    expect(validateWorkUnits(['folder' as any], [], undefined).length, 1);
    expect(validateWorkUnits(['series'], [{code: 'f', kind: 'study', connection: 'c', sql: 's'}], undefined)
      .some((e) => e.includes('does not enable')), true);
    expect(validateWorkUnits(['campaign'], [], ['hit']).some((e) => e.includes('needs the series')), true);
    expect(validateWorkUnits(['series'], [], ['hit', 'Hit']).length, 1, 'maturity values are case-blind');
  });

  test('units land under their program by code; unknown programs and duplicates are rejected', async () => {
    const p = plan({rows: [
      {program: 'P1', code: 'S-A', name: 'Amides', maturity: 'lead', status: 'active'},
      {program: 'P9', code: 'S-B', name: 'B'},
      {program: 'P1', code: 'S-A', name: 'again'},
      {program: 'P2', code: null, name: 'no code'},
    ], seriesMaturity: ['hit', 'lead']});
    expect(p.inserts.map((u) => `${u.program_id}:${u.code}:${u.maturity}`).join(), 'id-1:S-A:lead');
    expect(p.rejected.map((r) => r.row).join(), '2,3,4');
    expect(p.rejected[0].reason.includes('P9 is not a program of X'), true);
  });

  test('choices are checked, unknown columns kept as details, vocabulary drift warned', async () => {
    const p = plan({rows: [{program: 'P1', code: 'S-A', name: 'A', status: 'parked', maturity: 'mature', chemist: 'kim'}],
      seriesMaturity: ['hit', 'lead']});
    const [u] = p.inserts;
    expect(u.status === undefined, true, 'an invalid choice is left unset');
    expect(JSON.stringify(u.details), '{"chemist":"kim"}');
    expect(p.warnings.length, 2, p.warnings.join('; '));
  });

  test('a change is an update; the same values are unchanged; dates compare by instant', async () => {
    const existing = [unit('id-1', 'C1', {type: 'screening', status: 'active', started_on: '2026-01-12T00:00:00.000Z'})];
    const same = planWorkUnits({portfolio: 'X', kind: 'campaign', feed: 'f', programs, existing,
      rows: [{program: 'P1', code: 'C1', name: 'C1', type: 'screening', status: 'active', started_on: '2026-01-12'}]});
    expect(same.unchanged.join(), 'P1/C1', JSON.stringify(same.updates));
    const moved = planWorkUnits({portfolio: 'X', kind: 'campaign', feed: 'f', programs, existing,
      rows: [{program: 'P1', code: 'C1', name: 'C1', type: 'screening', status: 'closed'}]});
    expect(JSON.stringify(moved.updates[0].changes), '{"status":"closed"}');
  });

  test('a unit gone from its feed is flagged, not deleted; it comes back when the feed does', async () => {
    const gone = plan({existing: [unit('id-1', 'S-A', {}), unit('id-1', 'S-M', {}, {source: 'manual'})]});
    expect(gone.missing.map((m) => m.key).join(), 'P1/S-A', 'a manual unit is never missing');
    const back = plan({existing: [unit('id-1', 'S-A', {}, {missingFromFeed: true})],
      rows: [{program: 'P1', code: 'S-A', name: 'S-A'}]});
    expect(back.returned.map((r) => r.key).join(), 'P1/S-A');
    expect(back.missing.length, 0);
  });
});

let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  for (const code of ['EX-SMD', 'EX-BIOD', 'EX-BIOC']) {
    await syncPortfolio(code);
    await syncWorkUnits(code);
  }
})();

const unitsOf = async (portfolio: string, program: string, kind: any) => {
  const pf = await portfolioByCode(portfolio);
  const p = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${program}"`}))[0];
  return (await workUnitsOf(kind, [p.id])).sort((a, b) => a.code.localeCompare(b.code));
};

category('Portal: work units on the example portfolios', () => {
  before(async () => { await setUp(); });

  test('each portfolio enables its own units', async () => {
    expect(enabledWorkUnits(await portfolioByCode('EX-SMD')).join(), 'campaign,series');
    expect(enabledWorkUnits(await portfolioByCode('EX-BIOC')).join(), 'study,process_stage');
    expect(enabledWorkUnits(await portfolioByCode('EX-BIOD')).join(), 'campaign');
  });

  test('discovery: series and campaigns under their projects', async () => {
    const series = await unitsOf('EX-SMD', 'EXS-KIN1', 'series');
    expect(series.map((s) => `${s.code}:${s.maturity}:${s.status}`).join(), 'KIN-A:validated:active,KIN-B:hit:backup');
    const campaigns = await unitsOf('EX-SMD', 'EXS-KIN1', 'campaign');
    expect(campaigns.map((c) => `${c.code}:${c.type}:${c.outcome}`).join(),
      'DMTA-K1-07:dmta_cycle:open,HTS-K1:screening:succeeded', 'a campaign succeeds as a state');
  });

  test('development: studies complete, process stages in order', async () => {
    const studies = await unitsOf('EX-BIOC', 'EXB-201', 'study');
    expect(studies.map((s) => `${s.code}:${s.status}`).join(), 'STAB-201:running,TOX-201:completed');
    const stages = (await unitsOf('EX-BIOC', 'EXB-201', 'process_stage'))
      .sort((a, b) => a.process_type.localeCompare(b.process_type) || a.ordinal - b.ordinal);
    expect(stages.map((s) => s.code).join(), 'FORM,FILL,USP,DSP');
  });

  test('a campaign for a project the registry does not list is rejected and reported', async () => {
    const plans = await syncWorkUnits('EX-SMD', true);
    const campaigns = plans.find((p) => p.kind === 'campaign')!;
    expect(campaigns.rejected.length, 1);
    expect(campaigns.rejected[0].reason.includes('EXS-XYZ9'), true);
  });

  test('a second work-unit sync changes nothing', async () => {
    for (const code of ['EX-SMD', 'EX-BIOD', 'EX-BIOC']) {
      for (const p of await syncWorkUnits(code, true)) {
        expect(p.inserts.length + p.updates.length + p.missing.length + p.returned.length, 0,
          `${code} ${p.kind}: ${JSON.stringify(p.updates)}`);
      }
    }
  });

  test('a unit added by hand is never flagged by a sync', async () => {
    const pf = await portfolioByCode('EX-SMD');
    const kin = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "EXS-KIN1"`}))[0];
    const code = `TEST-WU-${Date.now() % 1e9}`;
    const [ins] = await portalDb.series.insert({program_id: kin.id, code, name: code});
    try {
      const got = await portalDb.series.get(ins.id);
      expect(got.source, 'manual');
      // a manual unit survives a sync of its kind
      await syncWorkUnits('EX-SMD');
      expect((await portalDb.series.get(ins.id)).missing_from_feed ?? false, false);
    } finally {
      await portalDb.series.delete(ins.id);
    }
  });
});
