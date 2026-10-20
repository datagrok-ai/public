import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext, readSources, recordDecision, routingAudit,
  readPortfolioFiles, setCriterion, shareProgram, stagesOf, syncPortfolio, verdictFor} from '../portal';
import {_package} from '../package-test';
import {withRestrictedUser} from './restricted-user';
import {buildPortfolioWorkspace} from '../workspace';
import {isParkingFolder, isRecordFolder, programFolderName} from '../model/workspace-plan';

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
    // a development gate on a research project must not place it: only the portal's own decisions do
    const disc = await programsOf('EX-BIOD');
    const log = await portalDb.stageTransitions.query({filter: disc.map((p) => `program_id = "${p.id}"`).join(' or '), limit: 10000});
    expect(log.filter((t) => t.decision === 'sync').length, 0, 'no registry transition in EX-BIOD');
    for (const p of disc) {
      const moves = log.filter((t) => t.program_id === p.id && ['advance', 'recycle'].includes(t.decision))
        .sort((a, b) => +new Date(a.decided_on as any) - +new Date(b.decided_on as any));
      expect(p.stage_id ?? null, moves.length ? moves[moves.length - 1].to_stage_id ?? null : null, `${p.code}: placed by its decisions`);
    }
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
    // program rows reach the portfolio's groups through their workspace folders
    await buildPortfolioWorkspace('EX-SMD');
    const smCount = (await programsOf('EX-SMD')).length;
    const ran = await withRestrictedUser('pportal', async (u) => {
      // the table is readable by everyone (schema grant); its rows are not
      expect(await u.asUser(() => portalDb.portfolios.count({})) >= 3, true, 'portfolios are visible to all');
      expect(await u.asUser(() => portalDb.programs.count({})), 0, 'no program before joining a portfolio');

      const viewers = await grok.dapi.groups.find(sm.viewers!);
      const personal = await grok.dapi.groups.find(u.group);
      await grok.dapi.groups.includeTo(personal, viewers);
      try {
        const all = await u.asUser(() => portalDb.programs.query({}));
        const seen = all.filter((p) => p.source === 'feed');
        expect(seen.length, smCount, `saw: ${seen.map((p) => p.code).join()}`);
        expect(all.every((p) => p.portfolio_id === sm.id), true, 'a small-molecule viewer saw another portfolio');
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

category('Portal: workspace', () => {
  before(async () => { await setUp(); });

  const ROOT = 'Example Small Molecule Discovery';
  const WORKSTREAMS = ['Chemistry', 'Biology', 'DMPK', 'CompChem'];
  const subspaces = async (id: string) =>
    await grok.dapi.spaces.id(id).children.ofTypes('Project', false).list({pageSize: 1000});

  test('a folder per program, its workstreams beneath; a second build creates nothing', async () => {
    const first = await buildPortfolioWorkspace('EX-SMD');
    expect(first.failed, 0, first.entries.filter((e) => e.status !== 'ok').map((e) => e.status).join('; '));
    const programs = await programsOf('EX-SMD');
    for (const p of programs) {
      expect(!!p.space_id, true, `${p.code} has no folder recorded`);
      const names = (await subspaces(p.space_id!)).map((s) => s.friendlyName ?? s.name).sort().join();
      expect(names, [...WORKSTREAMS].sort().join(), `${p.code}'s workstreams`);
    }
    const again = await buildPortfolioWorkspace('EX-SMD');
    expect(again.created, 0, again.entries.filter((e) => e.action === 'create').map((e) => e.path).join());
    // and the folders of the stages each program has reached (the example history places some)
    const smd = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-SMD')!;
    const stageCodes = new Map((await stagesOf(programs[0].stage_model_id)).map((s) => [s.id, s.ordinal]));
    const stageFolders = programs.map((p) => {
      const reached = smd.stageModel.stages.filter((s) => s.ordinal <= (stageCodes.get(p.stage_id ?? '') ?? 0));
      const paths = new Set(reached.flatMap((s) => s.folders ?? []).flatMap((f) => f.split('/').map((_, i, a) => a.slice(0, i + 1).join('/'))));
      return [...paths].filter((x) => !WORKSTREAMS.includes(x)).length;
    }).reduce((a, b) => a + b, 0);
    expect(again.reused, 1 + programs.length * (1 + WORKSTREAMS.length) + stageFolders);
    expect(again.notes.some((n) => n.startsWith('attributes written')), false, again.notes.join('; '));
  });

  test('a preview names what it would create and creates nothing', async () => {
    await buildPortfolioWorkspace('EX-SMD');
    const sm = await portfolioByCode('EX-SMD');
    const code = `WS-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: 'Preview only',
      stage_model_id: sm.stage_model_id, source: 'manual'});
    try {
      const preview = await buildPortfolioWorkspace('EX-SMD', true);
      const folder = `${ROOT}/${programFolderName({code, name: 'Preview only'})}`;
      expect(preview.entries.filter((e) => e.action === 'create').map((e) => e.path).sort().join(),
        [folder, ...WORKSTREAMS.map((w) => `${folder}/${w}`)].sort().join());
      expect((await portalDb.programs.get(ins.id)).space_id ?? null, null, 'a preview recorded a folder');
      const again = await buildPortfolioWorkspace('EX-SMD', true);
      expect(again.created, 1 + WORKSTREAMS.length, 'the first preview created something');
    } finally {
      await portalDb.programs.delete(ins.id);
    }
  });

  test('a program renamed in its registry keeps its folder', async () => {
    await buildPortfolioWorkspace('EX-SMD');
    const p = await programOf('EX-SMD', 'EXS-KIN1');
    await portalDb.programs.update(p.id, {name: `${p.name} (renamed)`});
    try {
      const r = await buildPortfolioWorkspace('EX-SMD');
      expect(r.created, 0, 'the rename made a second folder');
      expect((await portalDb.programs.get(p.id)).space_id, p.space_id);
      expect(r.notes.some((n) => n.startsWith('EXS-KIN1:') && n.includes('not renamed')), true, r.notes.join('; '));
    } finally {
      await portalDb.programs.update(p.id, {name: p.name});
    }
  });

  test('two builds under one root at once: the second is refused', async () => {
    const results = await Promise.allSettled([buildPortfolioWorkspace('EX-SMD'), buildPortfolioWorkspace('EX-SMD')]);
    expect(results.filter((r) => r.status === 'rejected').length, 1);
    expect(String((results.find((r) => r.status === 'rejected') as PromiseRejectedResult).reason)
      .includes('already running'), true);
  });

  test('the workspace is visible to the portfolio\'s viewers and to nobody else', async () => {
    await buildPortfolioWorkspace('EX-SMD');
    const sm = await portfolioByCode('EX-SMD');
    const ran = await withRestrictedUser('pportal', async (u) => {
      const roots = () => u.asUser(() => grok.dapi.spaces.filter(`name = "${ROOT}"`).list());
      expect((await roots()).length, 0, 'an outsider sees the workspace');
      const viewers = await grok.dapi.groups.find(sm.viewers!);
      const personal = await grok.dapi.groups.find(u.group);
      await grok.dapi.groups.includeTo(personal, viewers);
      try {
        const seen = await roots();
        expect(seen.length, 1, 'a viewer does not see the workspace');
        const all = await u.asUser(() => subspaces(seen[0].id));
        const ref = (f: any) => ({id: f.id, name: f.name, friendlyName: f.friendlyName ?? null});
        const folders = all.filter((f) => !isRecordFolder(ref(f)) && !isParkingFolder(ref(f)));
        expect(folders.length, (await programsOf('EX-SMD')).length, 'a viewer does not see the program folders');
        expect(all.some((f) => isRecordFolder(ref(f))), true, 'a viewer does not see _Sync history: it inherits the root\'s grants');
      } finally {
        await grok.dapi.groups.excludeFrom(personal, viewers);
      }
    });
    expect(ran !== null, true, 'no restricted session could be created on this stand');
  });
});
