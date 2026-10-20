import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {PortfolioConfig, validatePortfolio, validatePortfolioSet} from '../model/config';
import {
  diffCount, diffPortfolios, makeBundle, normalizePortfolio, planImport, readBundle, summarizeDiff, summarizeImport,
} from '../model/profile';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {completeCards, exportPortfolioConfig, importProfileBundle} from '../profile';
import {_package} from '../package-test';

const base = (): PortfolioConfig => ({
  code: 'P', name: 'P', setting: 'discovery', stageAuthority: 'portal', workspaceRoot: 'Root',
  workstreams: ['Chemistry', 'Biology'], workUnits: ['campaign'], schedule: 'every 6h',
  stageModel: {code: 'M', name: 'M', setting: 'discovery', stages: [
    {code: 'A', name: 'Alpha', ordinal: 1, isGate: false},
    {code: 'B', name: 'Beta', ordinal: 2, folders: ['DMPK'], criteria: [{code: 'c1', name: 'One', kind: 'document'}]},
  ]},
  catalogue: {unitOperations: [{code: 'OP1', name: 'Op one'}], models: [{code: 'Pkg:M1', name: 'M1'}]},
  feeds: [{code: 'f', system: 'S', connection: 'C', sql: 'select 1 as code, \'x\' as name'}],
});

const lines = (a: PortfolioConfig, b: PortfolioConfig) => {
  const d = diffPortfolios(a, b);
  return [...d.vocabulary, ...d.workspace, ...d.fields];
};

category('Portal: profile rules', () => {
  test('a default written out is no difference; empty lists are none', async () => {
    const a = base();
    const b: any = JSON.parse(JSON.stringify(a));
    b.programNoun = 'Program';
    b.stageModel.stages[1].isGate = true;
    b.stageModel.stages[1].criteria[0].required = true;
    b.catalogue.models[0].status = 'active';
    b.catalogue.models[0].package = 'Pkg';
    b.owners = [];
    expect(lines(a, b).join('; '), '');
    expect(JSON.stringify(normalizePortfolio(a)), JSON.stringify(normalizePortfolio(b)), 'the normal forms agree');
  });

  test('vocabulary: stages, criteria, work units and catalogue codes, in readable lines', async () => {
    const a = base();
    const b: any = JSON.parse(JSON.stringify(a));
    b.stageModel.stages.push({code: 'C', name: 'Gamma', ordinal: 3, criteria: [{code: 'c9', name: 'Nine', kind: 'decision'}]});
    b.stageModel.stages[1].name = 'Beta two';
    b.stageModel.stages[1].criteria[0].required = false;
    b.stageModel.stages[1].criteria.push({code: 'c2', name: 'Two', kind: 'decision'});
    b.workUnits = ['campaign', 'series'];
    b.catalogue.unitOperations = [{code: 'OP2', name: 'Op two'}];
    b.catalogue.models[0].card = {version: '1'};
    const d = diffPortfolios(a, b);
    expect(d.vocabulary.join(' | '), [
      'stage C (Gamma): added at ordinal 3 with criteria c9',
      'stage B.name: "Beta" → "Beta two"',
      'criterion B/c2: added',
      'criterion B/c1.required: true → false',
      'work units: added series',
      'unit operation OP1: removed',
      'unit operation OP2: added',
      'model Pkg:M1.card: set to {"approval":"draft","version":"1"}',
    ].join(' | '));
    expect(d.workspace.length + d.fields.length, 0);
  });

  test('workspace tree and other fields', async () => {
    const a = base();
    const b: any = JSON.parse(JSON.stringify(a));
    b.workstreams = ['Biology', 'DMPK'];
    b.stageModel.stages[1].folders = ['DMPK', 'DMPK/In vivo PK'];
    b.workspaceRoot = 'Root 2';
    b.schedule = '0 6 * * *';
    b.feeds[0].sql = 'select code, name from a_much_longer_registry_query_than_before where active = true order by code';
    b.groups = {viewers: 'V'};
    const d = diffPortfolios(a, b);
    expect(d.workspace.join(' | '), 'stage B folders: added DMPK/In vivo PK | workspace root: "Root" → "Root 2" | ' +
      'workstreams: added DMPK; removed Chemistry');
    expect(d.fields.join(' | '), 'feeds f.sql: changed (29 → 97 characters) | schedule: "every 6h" → "0 6 * * *" | ' +
      'groups: set to {"viewers":"V"}');
    expect(summarizeDiff('P', d).split('\n')[0], 'P: 6 difference(s)');
    expect(summarizeDiff('P', diffPortfolios(a, a)), 'P: no differences');
  });

  test('a bundle is checked as a set; import is refused when it does not validate', async () => {
    expect(readBundle('{"format": "x"}').errors.length > 0, true);
    expect(readBundle('not json').errors[0].startsWith('not JSON'), true);
    const a = base();
    const bundle = makeBundle([a], '2026-10-02T00:00:00Z');
    expect(readBundle(JSON.stringify(bundle)).errors.join('; '), '');
    expect(summarizeImport(planImport(bundle, new Map(), []), true).split('\n')[0],
      'Would import 1 portfolio(s): 1 new, 0 changed, 0 unchanged');
    const changed = {...a, description: 'changed'};
    const plan = planImport(makeBundle([changed], ''), new Map([['P', a]]), ['P']);
    expect(plan.portfolios[0].status, 'changed');
    const bad = makeBundle([{...a, stageAuthority: 'nobody' as any}], '');
    expect(summarizeImport(planImport(bad, new Map(), []), true).startsWith('Refused: the bundle does not validate'), true);
    // a handoff may go to a portfolio the instance already has
    const handing = {...a, stageModel: {...a.stageModel, stages: [...a.stageModel.stages.slice(0, 1),
      {...a.stageModel.stages[1], isTerminal: true}]}, handoff: {to: 'DEV', fromStage: 'B'}};
    expect(validatePortfolioSet([handing]).some((e) => e.includes('"DEV", which is not configured')), true);
    expect(validatePortfolioSet([handing], ['DEV']).join('; '), '');
  });
});

category('Portal: profile on the examples', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('export of each loaded example equals its file', async () => {
    const files = (await readPortfolioFiles(_package)).filter((c) => c.code.startsWith('EX-'));
    expect(files.map((f) => f.code).sort().join(), 'EX-BIOC,EX-BIOD,EX-SMD');
    for (const file of files) {
      const exported = await exportPortfolioConfig(file.code);
      expect(validatePortfolio(exported).join('; '), '', `${file.code} exported validates`);
      const d = diffPortfolios(exported, completeCards(file));
      expect(summarizeDiff(file.code, d), `${file.code}: no differences`);
    }
  });

  test('the functions: export, then diff a loaded portfolio against an edited file', async () => {
    const text: string = await grok.functions.call('ProjectPortal:exportPortfolio', {code: 'EX-SMD'});
    const edited = JSON.parse(text);
    expect(edited.code, 'EX-SMD');
    edited.workstreams = edited.workstreams.filter((w: string) => w !== 'CompChem');
    edited.stageModel.stages.find((s: any) => s.code === 'H2L').folders = ['Chemistry/SAR'];
    edited.schedule = 'every 12h';
    const out: string = await grok.functions.call('ProjectPortal:diffPortfolios', {a: 'EX-SMD', b: JSON.stringify(edited)});
    expect(out.split('\n')[0], 'EX-SMD (loaded) → EX-SMD (file): 3 difference(s)');
    expect(out.includes('workstreams: removed CompChem') && out.includes('stage H2L folders: added Chemistry/SAR') &&
      out.includes('schedule: "every 6h" → "every 12h"'), true, out);
  });

  test('a profile moves the examples: dry run, refusal, then an import that changes one field', async () => {
    const bundle: string = await grok.functions.call('ProjectPortal:exportProfile', {codes: 'EX-SMD,EX-BIOD,EX-BIOC'});
    const dry: string = await grok.functions.call('ProjectPortal:importProfile', {bundle, dryRun: true});
    expect(dry.split('\n')[0], 'Would import 3 portfolio(s): 0 new, 0 changed, 3 unchanged');
    const b = JSON.parse(bundle);
    const bad = JSON.parse(bundle);
    bad.portfolios[0].stageModel.stages[0].code = bad.portfolios[0].stageModel.stages[1].code;
    expect((await importProfileBundle(bad, false)).report.startsWith('Refused: the bundle does not validate'), true);
    const biod = b.portfolios.find((p: any) => p.code === 'EX-BIOD');
    const was = biod.description;
    biod.description = `${was} (imported)`;
    const preview = await importProfileBundle(b, true);
    expect(preview.plan.portfolios.find((p) => p.code === 'EX-BIOD')!.status, 'changed');
    expect((await portfolioByCode('EX-BIOD')).description, was, 'a dry run writes nothing');
    try {
      const r = await importProfileBundle(b, false);
      expect(r.report.startsWith('Imported 3 portfolio(s): 0 new, 1 changed, 2 unchanged'), true, r.report);
      expect((await portfolioByCode('EX-BIOD')).description, `${was} (imported)`);
    } finally {
      await loadAllPortfolios(_package);
    }
    expect(diffCount(diffPortfolios(await exportPortfolioConfig('EX-BIOD'),
      (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-BIOD')!)), 0, 'the files put back');
  });
});
