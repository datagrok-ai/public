import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext, readSources, recordDecision, routingAudit,
  setCriterion, shareProgram, stagesOf, syncPortfolio, verdictFor} from '../portal';
import {_package} from '../package-test';
import {withRestrictedUser} from './restricted-user';

/** The three example portfolios, loaded and synced once per run. Their feeds are synthetic
 * registries written as VALUES on System:Datagrok, so every stand has them. Discovery syncs before
 * development so lineage resolves in one pass. */
const PORTFOLIOS = ['EX-SMD', 'EX-BIOD', 'EX-BIOC'];
let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  for (const code of PORTFOLIOS)
    await syncPortfolio(code);
})();

const programsOf = async (code: string) =>
  portalDb.programs.query({filter: `portfolio_id = "${(await portfolioByCode(code)).id}" and source = "feed"`});
const programOf = async (portfolio: string, code: string) =>
  (await programsOf(portfolio)).find((p) => p.code === code)!;

category('Portal: portfolios on one instance', () => {
  before(async () => { await setUp(); });

  test('three portfolios, three stage models', async () => {
    const models = new Set<string>();
    for (const code of PORTFOLIOS)
      models.add((await portfolioByCode(code)).stage_model_id);
    expect(models.size, 3);
    const dev = await stagesOf((await portfolioByCode('EX-BIOC')).stage_model_id);
    expect(dev.map((s) => s.code).join(), 'D1,D2,D3,D4,D5,D6,D7');
  });

  test('the same stage code lives in two models without colliding', async () => {
    const lo = await portalDb.stages.query({filter: 'code = "LO"'});
    const models = new Set(lo.map((s) => s.stage_model_id));
    expect(models.size >= 2, true, `LO rows: ${lo.length}`);
  });

  test('loading is idempotent', async () => {
    const before = await portalDb.stages.count({});
    await loadAllPortfolios(_package);
    expect(await portalDb.stages.count({}), before);
  });

  test('each feed takes its own slice of a shared registry', async () => {
    const codes = async (p: string) => (await programsOf(p)).map((x) => x.code).sort().join();
    expect(await codes('EX-SMD'), 'EXS-EPI5,EXS-GPC3,EXS-ION4,EXS-KIN1,EXS-PRO2');
    expect(await codes('EX-BIOD'), 'EXB-101,EXB-102,EXB-103,EXB-104,EXB-105');
    expect(await codes('EX-BIOC'), 'EXB-104,EXB-105,EXB-201,EXB-202,EXB-203,EXB-204,EXB-205');
  });

  test('stages follow the authority: the registry places development, the portal places discovery', async () => {
    const dev = await programsOf('EX-BIOC');
    const stageOf = async (code: string) => {
      const p = dev.find((x) => x.code === code)!;
      return p.stage_id ? (await portalDb.stages.get(p.stage_id)).code : null;
    };
    expect(await stageOf('EXB-201'), 'D4', 'the alias GD4 resolves');
    expect(await stageOf('EXB-104'), 'D1', 'the alias Gate D1 resolves');
    expect(await stageOf('EXB-205'), null, 'an unresolvable spelling leaves the program unstaged');
    const disc = await programsOf('EX-BIOD');
    expect(disc.filter((p) => p.stage_id).length, 0, 'a development gate on a research project must not place it');
    expect(disc.find((p) => p.code === 'EXB-101')?.source_stage, '3');
  });

  test('a molecule handed on is a new program linked to its predecessor', async () => {
    const research = await programOf('EX-BIOD', 'EXB-104');
    const dev = await programOf('EX-BIOC', 'EXB-104');
    expect(dev.predecessor_id, research.id);
    expect(research.status, 'transitioned', 'nominated in the registry, handed on in the portal');
  });

  test('a molecule still active in research while in development is an overlap', async () => {
    const research = await programOf('EX-BIOD', 'EXB-105');
    const dev = await programOf('EX-BIOC', 'EXB-105');
    expect(dev.predecessor_id, research.id);
    expect(research.status, 'active');
  });

  test('an unresolvable stage spelling is reported on every sync', async () => {
    const [plan] = await syncPortfolio('EX-BIOC', true);
    expect(plan.unresolvedStages.map((u) => `${u.code}=${u.spelling}`).join(), 'EXB-205=Gate 9z');
  });

  test('a second sync changes nothing', async () => {
    for (const code of PORTFOLIOS) {
      const [plan] = await syncPortfolio(code, true);
      expect(plan.inserts.length + plan.updates.length + plan.missing.length, 0,
        `${code}: ${JSON.stringify({i: plan.inserts.map((x) => x.code), u: plan.updates})}`);
    }
  });

  test('the routing audit finds what no portfolio claims, and what two claim', async () => {
    const reports = await routingAudit(await readSources(_package));
    const bio = reports.find((r) => r.system === 'example-bio-registry')!;
    expect(bio.unclaimed.join(), 'EXB-301,EXB-302', 'a development small molecule and an unclassified project');
    expect(bio.claimedTwice.map((c) => c.code).sort().join(), 'EXB-104,EXB-105');
    expect(bio.total, 12);
    const sm = reports.find((r) => r.system === 'example-sm-registry')!;
    expect(sm.unclaimed.length, 0, sm.unclaimed.join());
  });
});

category('Portal: gate decisions', () => {
  before(async () => { await setUp(); });

  test('a portal-owned gate: criteria, advance, log', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const code = `TEST-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, source: 'manual'});
    try {
      // unstaged → TV: no criteria on entry
      await recordDecision(ins.id, 'advance', null, 'kick-off');
      let ctx = await programContext(ins.id);
      const tv = ctx.stages.find((s) => s.code === 'TV')!;
      expect(ctx.program.stage_id, tv.id);
      // TV → HF needs the validated assay
      let v = verdictFor(ctx, 'advance');
      expect(v.allowed, false);
      expect(v.outstanding.join(), 'assay_ready');
      const hf = ctx.stages.find((s) => s.code === 'HF')!;
      const assay = ctx.criteria.find((c) => c.stageId === hf.id && c.code === 'assay_ready')!;
      await setCriterion(ins.id, assay.id, 'met', 'Z\' = 0.71');
      ctx = await programContext(ins.id);
      v = verdictFor(ctx, 'advance');
      expect(v.allowed, true, v.reasons.join('; '));
      await recordDecision(ins.id, 'advance', null, 'assay validated');
      const log = await portalDb.stageTransitions.query({filter: `program_id = "${ins.id}"`});
      expect(log.length, 2);
      expect(log.every((t) => t.decision === 'advance' && t.decided_by === grok.shell.user.id), true);
      // a stage of the development model is not reachable from a discovery program
      const devStage = (await stagesOf((await portfolioByCode('EX-BIOC')).stage_model_id))[0];
      expect(verdictFor(await programContext(ins.id), 'advance', devStage.id).allowed, false);
    } finally {
      await portalDb.programs.delete(ins.id);
    }
  });

  test('a registry-owned gate cannot be moved in the portal', async () => {
    const alpha = await programOf('EX-BIOC', 'EXB-201');
    const v = verdictFor(await programContext(alpha.id), 'advance');
    expect(v.allowed, false);
    expect(v.reasons.some((r) => r.includes('registry owns')), true);
    let refused = false;
    try {
      await recordDecision(alpha.id, 'advance', null, 'try');
    } catch {
      refused = true;
    }
    expect(refused, true);
  });

  test('a feed stage change is logged as a sync transition', async () => {
    const alpha = await programOf('EX-BIOC', 'EXB-201');
    const log = await portalDb.stageTransitions.query({filter: `program_id = "${alpha.id}"`});
    expect(log.some((t) => t.decision === 'sync'), true);
  });
});

category('Portal: isolation between portfolios', () => {
  before(async () => { await setUp(); });

  test('a user sees only the portfolios whose groups they are in', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const smCount = (await programsOf('EX-SMD')).length;
    const ran = await withRestrictedUser('pportal', async (u) => {
      // the table is readable by everyone (schema grant); its rows are not
      expect(await u.asUser(() => portalDb.portfolios.count({})) >= 3, true, 'portfolios are visible to all');
      expect(await u.asUser(() => portalDb.programs.count({})), 0, 'no program before joining a portfolio');

      const viewers = await grok.dapi.groups.find(sm.viewers!);
      const personal = await grok.dapi.groups.find(u.group);
      await grok.dapi.groups.includeTo(personal, viewers);
      try {
        const seen = await u.asUser(() => portalDb.programs.query({}));
        expect(seen.length, smCount, `saw: ${seen.map((p) => p.code).join()}`);
        expect(seen.every((p) => p.portfolio_id === sm.id), true, 'a small-molecule viewer saw another portfolio');
      } finally {
        await grok.dapi.groups.excludeFrom(personal, viewers);
      }
    });
    expect(ran !== null, true, 'no restricted session could be created on this stand');
  });

  test('a manually created program stays private until shared', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const code = `PRIV-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, source: 'manual'});
    try {
      await withRestrictedUser('pportal', async (u) => {
        const viewers = await grok.dapi.groups.find(sm.viewers!);
        const personal = await grok.dapi.groups.find(u.group);
        await grok.dapi.groups.includeTo(personal, viewers);
        try {
          const find = () => u.asUser(() => portalDb.programs.query({filter: `code = "${code}"`}));
          expect((await find()).length, 0, 'an unshared program leaked to the portfolio');
          await shareProgram(ins.id, sm);
          expect((await find()).length, 1, 'sharing with the portfolio did not reach its viewers');
        } finally {
          await grok.dapi.groups.excludeFrom(personal, viewers);
        }
      });
    } finally {
      await portalDb.programs.delete(ins.id);
    }
  });
});
