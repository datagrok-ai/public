import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {checkDesignAccess, describeDesignSet, designSetCode, proposalColumns, proposalsFromRows, rankProposals,
  selectionCompounds, summarizeDesignSet, validateDesignSet, validateSelectionSet} from '../model/design-sets';
import {createSelectionSet, designSetsOf, loadExampleDesignSets, selectionSetsOf} from '../design-sets';
import {openDesignSetDialog} from '../views/design-set-dialog';
import {filterValue, ROW_LIMIT} from '../rows';
import {_package} from '../package-test';

category('Portal: design set rules', () => {
  test('a design set extends a series of the program; labels once, a structure each', async () => {
    const ok = {name: 'Para', series: 'KIN-A', source: 'manual' as const, proposals: [{label: 'D1', structure: 'CCO', score: 1}]};
    expect(validateDesignSet(ok, ['KIN-A']).join(), '');
    expect(validateDesignSet({...ok, series: 'KIN-Z'}, ['KIN-A']).join(), 'the program has no series KIN-Z');
    expect(validateDesignSet({...ok, source: 'ai' as any, proposals: []}, ['KIN-A']).join('; '),
      'source is one of enumerated, generative, manual; a design set proposes at least one structure');
    expect(validateDesignSet({...ok, proposals: [{label: 'D1', structure: 'C'}, {label: 'D1', structure: ''}]}, ['KIN-A']).join('; '),
      'label D1 is used twice; proposal D1 needs a structure');
    expect(checkDesignAccess({isApprover: false, isContributor: true}).length, 0);
    expect(checkDesignAccess({isApprover: true, isContributor: true, programStatus: 'terminated'}).join(), 'The program is terminated');
  });

  test('ranked by score, ties share a rank, unscored last; lower is better when asked', async () => {
    const p = [{label: 'a', score: 7}, {label: 'b', score: 8}, {label: 'c', score: 7}, {label: 'd', score: null}, {label: 'e', score: 6}];
    expect([...rankProposals(p).entries()].map(([k, v]) => `${k}${v ?? '-'}`).join(), 'b1,a2,c2,e4,d-');
    expect([...rankProposals(p, false).entries()].map(([k, v]) => `${k}${v ?? '-'}`).join(), 'e1,a2,c2,b4,d-');
    const s = summarizeDesignSet([{label: 'x', rank: 2, chosen: true}, {label: 'y', rank: 1}, {label: 'z'}]);
    expect(describeDesignSet(s), '3 proposed · 1 chosen for synthesis · best y');
  });

  test('proposals from a table: the structure, label and score columns; labels made up when missing', async () => {
    const cols = proposalColumns([{name: 'ID'}, {name: 'Mol', semType: 'Molecule'}, {name: 'Score', numeric: true}]);
    expect(`${cols.structure} ${cols.label} ${cols.score}`, 'Mol ID Score');
    expect(proposalColumns([{name: 'smiles'}, {name: 'pIC50'}], 'pIC50').score, 'pIC50');
    const p = proposalsFromRows([{smiles: 'CCO', s: '7.1'}, {smiles: 'CCN', s: ''}], {structure: 'smiles', label: null, score: 's'}, 'KIN-A-D');
    expect(p.map((x) => `${x.label}:${x.structure}:${x.score}`).join(), 'KIN-A-D-01:CCO:7.1,KIN-A-D-02:CCN:null');
    expect(designSetCode('KIN-A', new Date('2026-10-03T05:00:00Z'), ['KIN-A-2026-10-03']), 'KIN-A-2026-10-03-2');
  });

  test('a selection set: a rationale, the program\'s own compounds, once each', async () => {
    const ok = {name: 'Hits', campaign: 'HTS-K1', rationale: 'potent', compounds: ['C1', 'C2']};
    expect(validateSelectionSet(ok, ['HTS-K1'], ['C1', 'C2']).join(), '');
    expect(validateSelectionSet({...ok, campaign: 'X', rationale: ''}, ['HTS-K1'], ['C1', 'C2']).join('; '),
      'the program has no campaign X; a selection set needs a rationale: why these compounds');
    expect(validateSelectionSet({...ok, compounds: ['C1', 'C9', 'C1']}, [], ['C1']).join('; '),
      'the program has no campaign HTS-K1; not compounds of the program: C9; a compound is listed twice');
    expect(selectionCompounds({items: ['C1', ' ', 'C2']}).join(), 'C1,C2');
  });
});

category('Portal: design sets', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  const undo: (() => Promise<any>)[] = [];

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`, limit: ROW_LIMIT}))[0];
    await loadExampleDesignSets(_package);
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
    document.querySelectorAll('[data-u2-name="design-set-dialog"]').forEach((d) => d.remove());
  });

  test('the examples on EXS-KIN1: two KIN-A sets, ranked, chosen ones marked; loaded once', async () => {
    const sets = await designSetsOf(kin1.id);
    const para = sets.find((s) => s.set.code === 'KIN-A-2026-04-20')!;
    expect(`${para.set.series_code} ${para.set.source} ${para.proposals.length}`, 'KIN-A enumerated 6');
    expect(para.proposals.map((p) => `${p.label}:${p.rank}`).slice(0, 2).join(), 'KIN-A-D01:1,KIN-A-D02:2');
    expect(describeDesignSet(summarizeDesignSet(para.proposals)), '6 proposed · 4 chosen for synthesis · best KIN-A-D01');
    expect(sets.some((s) => s.set.code === 'KIN-A-2026-09-01'), true);
    const again = await loadExampleDesignSets(_package);
    expect(again.includes('design set KIN-A-2026-04-20 already there'), true, again);
    expect((await selectionSetsOf(kin1.id)).find((x) => x.code === 'HTS-K1-HITS')?.campaign_code, 'HTS-K1');
  });

  test('a set from a table through the function; chosen and unchosen; the table', async () => {
    const df = DG.DataFrame.fromCsv('label,smiles,score\nT1,Nc1nccc(-c2ccc(I)cc2)n1,6.5\nT2,Nc1nccc(-c2ccc(O)cc2)n1,7.5');
    const name = `Test set ${Date.now()}`;
    const report: string = await grok.functions.call('ProjectPortalCore:createDesignSet', {portfolio: 'EX-SMD', program: 'EXS-KIN1',
      series: 'KIN-B', name, source: 'manual', proposals: df});
    const set = (await designSetsOf(kin1.id)).find((s) => s.set.name === name)!;
    undo.push(() => portalDb.designSets.delete(set.set.id));
    expect(report.endsWith('2 proposed · 0 chosen for synthesis · best T2'), true, report);
    const chose: string = await grok.functions.call('ProjectPortalCore:chooseForSynthesis', {portfolio: 'EX-SMD', program: 'EXS-KIN1',
      designSet: set.set.code, labels: 'T1, T2'});
    expect(chose.endsWith('2 proposal(s) chosen for synthesis'), true, chose);
    await grok.functions.call('ProjectPortalCore:chooseForSynthesis', {portfolio: 'EX-SMD', program: 'EXS-KIN1',
      designSet: set.set.code, labels: 'T1', chosen: false});
    const back = (await designSetsOf(kin1.id)).find((s) => s.set.id === set.set.id)!;
    expect(back.proposals.map((p) => `${p.label}:${p.chosen}:${p.chosen_by === DG.User.current().id}`).join(), 'T2:true:true,T1:false:true');
    const table: DG.DataFrame = await grok.functions.call('ProjectPortalCore:designSets', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
    expect(table.col('Structure')!.semType, DG.SEMTYPE.MOLECULE);
    expect(table.col('Label')!.toList().includes('T2'), true);
  });

  test('a selection set of another program\'s compound is refused', async () => {
    await expectExceptionAsync(async () => {
      await createSelectionSet(kin1, {name: 'x', rationale: 'y', compounds: ['EXC-2001']});
    });
  });

  test('the dialog: proposals by rank with their choice', async () => {
    const para = (await designSetsOf(kin1.id)).find((s) => s.set.code === 'KIN-A-2026-04-20')!;
    const dlg = openDesignSetDialog(kin1, para, []);
    try {
      expect(dlg.root.querySelectorAll('[data-u2-name^="proposal-"]').length, 6);
      expect((dlg.root.querySelector('[data-u2-name="choose-KIN-A-D01"] input') as HTMLInputElement).checked, true);
      expect((dlg.root.querySelector('[data-u2-name="choose-KIN-A-D06"] input') as HTMLInputElement).checked, false);
      expect(!!dlg.root.querySelector('[data-u2-name="save-choices-button"]'), true);
    } finally {
      dlg.close();
    }
  });
});
