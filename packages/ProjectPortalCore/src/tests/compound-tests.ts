import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {CompoundPlanInput, CompoundRef, LinkRef, planCompounds} from '../model/compounds';
import {portalDb} from '../generated/db';
import {compoundsOf, loadAllPortfolios, portfolioByCode, setCompoundRole, syncCompounds, syncPortfolio} from '../portal';
import {_package} from '../package-test';

const programs = new Map([['P1', 'p1'], ['P2', 'p2']]);
const compounds = (codes: string[]) => new Map<string, CompoundRef>(codes.map((c) =>
  [c, {id: `id-${c}`, code: c, values: {name: null, smiles: 'C', aliases: null, batches: null, details: null}}]));
const plan = (over: Partial<CompoundPlanInput>) => planCompounds({portfolio: 'X', system: 'reg', rows: [], programs,
  compounds: new Map(), links: [], ...over});
const link = (programId: string, code: string, role: any, extra: Partial<LinkRef> = {}): LinkRef =>
  ({id: `l-${programId}-${code}`, programId, compoundId: `id-${code}`, role, source: 'feed', missingFromFeed: false, ...extra});

category('Portal: compound rules', () => {
  test('a compound is one row however many programs hold it; each program gets its link and role', async () => {
    const p = plan({rows: [
      {program: 'P1', compound: 'C1', smiles: 'CCO', role: 'lead', batches: 'C1-01|C1-02|C1-01'},
      {program: 'P2', compound: 'C1', smiles: 'CCO', role: 'tool'},
    ]});
    expect(p.compoundInserts.length, 1);
    expect(p.compoundInserts[0].batches, 'C1-01|C1-02', 'batches are a set');
    expect(p.linkInserts.map((l) => `${l.programCode}:${l.role}`).join(), 'P1:lead,P2:tool');
  });

  test('dirty rows: unknown program, duplicate link, missing code rejected; unknown role taken as explored', async () => {
    const p = plan({rows: [{program: 'P9', compound: 'C1'}, {program: 'P1', compound: 'C1', role: 'Star'},
      {program: 'P1', compound: 'C1'}, {program: 'P1', compound: ''}]});
    expect(p.rejected.map((r) => r.row).join(), '1,3,4');
    expect(p.linkInserts[0].role, 'explored');
    expect(p.warnings.some((w) => w.includes('"Star"')), true);
  });

  test('a role the registry moves is a logged change; the same role is unchanged', async () => {
    const p = plan({compounds: compounds(['C1', 'C2']), links: [link('p1', 'C1', 'hit'), link('p1', 'C2', 'lead')],
      rows: [{program: 'P1', compound: 'C1', smiles: 'C', role: 'lead'}, {program: 'P1', compound: 'C2', smiles: 'C', role: 'lead'}]});
    expect(p.roleChanges.map((c) => `${c.compoundCode}:${c.from}>${c.to}`).join(), 'C1:hit>lead');
    expect(p.unchanged, 1);
    expect(p.compoundUpdates.length, 0);
  });

  test('a role set by hand is kept against the registry, and reported', async () => {
    const p = plan({compounds: compounds(['C1']), links: [link('p1', 'C1', 'candidate', {source: 'manual'})],
      rows: [{program: 'P1', compound: 'C1', smiles: 'C', role: 'backup'}]});
    expect(p.roleChanges.length, 0);
    expect(p.warnings.join().includes('kept candidate, set by hand'), true);
  });

  test('a link the feed drops is flagged and comes back; a manual link is never missing', async () => {
    const gone = plan({compounds: compounds(['C1', 'C2']),
      links: [link('p1', 'C1', 'lead'), link('p1', 'C2', 'tool', {source: 'manual'})]});
    expect(gone.linkMissing.map((m) => m.key).join(), 'P1/C1');
    const back = plan({compounds: compounds(['C1']), links: [link('p1', 'C1', 'lead', {missingFromFeed: true})],
      rows: [{program: 'P1', compound: 'C1', smiles: 'C', role: 'lead'}]});
    expect(back.linkReturned.length, 1);
  });
});

let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  await syncPortfolio('EX-SMD');
  await syncCompounds('EX-SMD');
})();

const program = async (code: string) => {
  const sm = await portfolioByCode('EX-SMD');
  return (await portalDb.programs.query({filter: `portfolio_id = "${sm.id}" and code = "${code}"`}))[0];
};

category('Portal: compounds on the example portfolios', () => {
  before(async () => { await setUp(); });

  test('a project holds its compounds with their roles', async () => {
    const kin = await compoundsOf((await program('EXS-KIN1')).id);
    const roles = (r: string) => kin.filter((x) => x.link.role === r).length;
    expect(kin.length, 12);
    expect(`${roles('lead')}/${roles('backup')}/${roles('explored')}/${roles('discontinued')}/${roles('tool')}`, '6/2/2/1/1');
    expect(kin.find((x) => x.compound.code === 'EXC-1001')!.compound.batches, 'EXC-1001-01|EXC-1001-02');
  });

  test('a tool compound serves two projects: one compound, two links, two roles', async () => {
    const t1 = await portalDb.compounds.query({filter: 'system = "example-sm-compounds" and code = "EXC-3001"'});
    expect(t1.length, 1);
    const links = await portalDb.programCompounds.query({filter: `compound_id = "${t1[0].id}"`});
    expect(links.map((l) => l.role).sort().join(), 'lead,tool');
  });

  test('each link starts its history; a compound for an unknown project is rejected', async () => {
    const kin = await program('EXS-KIN1');
    const log = await portalDb.compoundRoleChanges.query({filter: `program_id = "${kin.id}"`});
    expect(log.length >= 12, true, `${log.length} history rows`);
    const plan = (await syncCompounds('EX-SMD', true))!;
    expect(plan.rejected.map((r) => r.reason).join().includes('EXS-XYZ9'), true);
  });

  test('a second compound sync changes nothing', async () => {
    const plan = (await syncCompounds('EX-SMD', true))!;
    expect(plan.compoundInserts.length + plan.compoundUpdates.length + plan.linkInserts.length +
      plan.roleChanges.length + plan.linkMissing.length, 0, JSON.stringify(plan.compoundUpdates.slice(0, 2)));
  });

  test('a role changed by hand is logged with its rationale', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const code = `TEST-CPD-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, source: 'manual'});
    try {
      const c = (await portalDb.compounds.query({filter: 'system = "example-sm-compounds" and code = "EXC-3003"'}))[0];
      await setCompoundRole(ins.id, c.id, 'backup', 'profiled as a backup');
      await setCompoundRole(ins.id, c.id, 'candidate', 'cleanest safety profile');
      const log = (await portalDb.compoundRoleChanges.query({filter: `program_id = "${ins.id}"`}))
        .sort((a, b) => String(a.changed_on).localeCompare(String(b.changed_on)));
      expect(log.map((l) => `${l.from_role ?? '-'}>${l.to_role}`).join(), '->backup,backup>candidate');
      expect(log[1].rationale, 'cleanest safety profile');
      let refused = false;
      await setCompoundRole(ins.id, c.id, 'lead', ' ').catch(() => refused = true);
      expect(refused, true, 'a role change needs a rationale');
    } finally {
      await portalDb.programs.delete(ins.id);
    }
  });
});
