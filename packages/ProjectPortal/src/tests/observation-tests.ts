import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {checkObservationAccess, checkRetraction, describeObservation, ObservationView, observationsAbout, sortObservations,
  subjectOptions, validateObservation} from '../model/observations';
import {FeedRecords, filterFeed, programFeedEntries} from '../model/feed';
import {addObservation, observationsOf, observationTargetsOf, retractObservation} from '../observations';
import {programFeedOf} from '../feed';
import {openObservationDialog} from '../views/observation-dialog';
import {filterValue, ROW_LIMIT} from '../rows';
import {_package} from '../package-test';

const targets = {series: ['KIN-B', 'KIN-A'], compound: ['EXC-1001'], campaign: ['HTS-K1']};
const view = (id: string, over: Partial<ObservationView> = {}): ObservationView => ({id, subject: 'series', code: 'KIN-A',
  title: `T${id}`, finding: 'f', evidence: [], status: 'current', by: 'u1', on: '2026-10-01T00:00:00Z', ...over});

category('Portal: observation rules', () => {
  test('about the program or something it holds; a title, a finding, evidence outside the gate', async () => {
    expect(subjectOptions(targets).map((o) => o.label).join(), 'The program,series KIN-A,series KIN-B,compound EXC-1001,campaign HTS-K1');
    const ok = {subject: 'series' as const, code: 'KIN-A', title: 'Para wins', finding: 'Potency up, solubility flat'};
    expect(validateObservation(ok, targets).join(), '');
    expect(validateObservation({...ok, subject: 'program', code: null}, targets).join(), '');
    expect(validateObservation({...ok, code: 'KIN-Z'}, targets).join(), 'the program has no series KIN-Z');
    expect(validateObservation({...ok, subject: 'study', code: 'S1'}, targets).join(), 'the program has no study S1');
    expect(validateObservation({...ok, code: ''}, targets).join(), 'an observation on a series names it');
    expect(validateObservation({...ok, title: ' ', finding: ''}, targets).join('; '),
      'an observation needs a title; an observation needs the finding itself');
    expect(validateObservation({...ok, evidence: [{kind: 'criterion', target: 'c1', name: 'C'}]}, targets)[0],
      'evidence of an observation is one of model_run, snapshot, file, url, not criterion');
    expect(validateObservation({...ok, evidence: [{kind: 'url', target: 'nope', name: 'x'}]}, targets).join(),
      '"nope" is not an http(s) address');
  });

  test('who records and retracts; current first, newest first; narrowed to a subject', async () => {
    expect(checkObservationAccess({isApprover: false, isContributor: true}).length, 0);
    expect(checkObservationAccess({isApprover: false, isContributor: false}).length, 1);
    expect(checkRetraction({status: 'current', by: 'u1'}, 'u1', false, 'superseded').join(), '');
    expect(checkRetraction({status: 'current', by: 'u1'}, 'u2', false, 'x').join(), 'Only its author or an approver retracts an observation');
    expect(checkRetraction({status: 'current', by: 'u1'}, 'u2', true, ' ').join(), 'Retracting needs a reason');
    expect(checkRetraction({status: 'retracted', by: 'u1'}, 'u1', false, 'x').join(), 'It is retracted already');
    const list = [view('1', {on: '2026-09-01T00:00:00Z'}), view('2', {status: 'retracted', on: '2026-10-02T00:00:00Z'}),
      view('3', {on: '2026-09-15T00:00:00Z', subject: 'compound', code: 'EXC-1001'})];
    expect(sortObservations(list).map((o) => o.id).join(), '3,1,2');
    expect(observationsAbout(list, 'series', 'KIN-A').map((o) => o.id).join(), '1,2');
    expect(observationsAbout(list).length, 3);
    expect(describeObservation(list[1]), 'series KIN-A · T2 (retracted)');
  });

  test('in the feed: recorded and retracted, tagged with the series or the compound', async () => {
    const r: FeedRecords = {program: {id: 'p1', code: 'P1'}, stages: [], decisions: [], evidence: [], runs: [], selections: [],
      documents: [], snapshots: [], roleChanges: [], syncs: [], news: [], observations: [
        {id: 'o1', subject: 'series', code: 'KIN-A', title: 'Para wins', by: 'u1', on: '2026-10-01T10:00:00Z',
          retractedBy: 'u2', retractedOn: '2026-10-02T10:00:00Z', retraction: 'new data'},
        {id: 'o2', subject: 'compound', code: 'EXC-1001', title: 'Clean hERG', by: 'u1', on: '2026-10-01T11:00:00Z'},
        {id: 'o3', subject: 'program', code: null, title: 'Kinase selectivity is the risk', by: 'u1', on: '2026-10-01T12:00:00Z'}]};
    const e = programFeedEntries(r);
    expect(e.map((x) => x.text).sort().join(' | '), 'Observation on compound EXC-1001: Clean hERG | Observation on series KIN-A: Para wins | ' +
      'Observation retracted on series KIN-A: Para wins — new data | Observation: Kinase selectivity is the risk');
    expect(e.every((x) => x.kind === 'observation'), true);
    const series = filterFeed(e, {series: 'KIN-A', seriesCompounds: ['EXC-1001']});
    expect(series.length, 3, 'the series\' two entries and its compound\'s');
  });
});

category('Portal: observations', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  const made: string[] = [];

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(smd.id)} and code = "EXS-KIN1"`, limit: ROW_LIMIT}))[0];
  });

  after(async () => {
    for (const id of made.splice(0))
      await portalDb.observations.delete(id).catch(() => {});
    document.querySelectorAll('[data-u2-name="observation-dialog"]').forEach((d) => d.remove());
  });

  test('EXS-KIN1 offers its series, campaigns and compounds', async () => {
    const t = await observationTargetsOf(smd, kin1.id);
    expect(t.series!.slice().sort().join(), 'KIN-A,KIN-B');
    expect(t.campaign!.includes('DMTA-K1-07'), true);
    expect(t.compound!.includes('EXC-1001'), true);
    expect(t.study === undefined, true, 'EX-SMD does not enable studies');
  });

  test('recorded on series KIN-A with a link: listed, in the function\'s table and the series\' feed', async () => {
    const o = await addObservation(kin1.id, {subject: 'series', code: 'KIN-A', title: `Para wins ${Date.now()}`,
      finding: 'para-F and para-Cl hold potency; solubility follows lipophilicity',
      evidence: [{kind: 'url', target: 'https://example.org/kin-a-review', name: 'KIN-A review'}]});
    made.push(o.id);
    expect(`${o.status} ${o.by === DG.User.current().id} ${o.evidence[0]?.name}`, 'current true KIN-A review');
    expect((await observationsOf(kin1.id))[0].id, o.id, 'newest current first');
    const df: DG.DataFrame = await grok.functions.call('ProjectPortal:observations',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', subject: 'series', code: 'KIN-A'});
    expect(df.col('Title')!.toList().includes(o.title), true);
    expect(df.col('About')!.toList().every((x) => x === 'series KIN-A'), true);
    const feed = await programFeedOf(smd, kin1, {series: 'KIN-A'});
    expect(feed.some((x) => x.kind === 'observation' && x.text === `Observation on series KIN-A: ${o.title}`), true);
  });

  test('refused on something the program does not hold; retracted with a reason, once', async () => {
    await expectExceptionAsync(async () => {
      await addObservation(kin1.id, {subject: 'compound', code: 'EXC-2001', title: 'x', finding: 'y'});
    });
    const report: string = await grok.functions.call('ProjectPortal:addObservation', {portfolio: 'EX-SMD', program: 'EXS-KIN1',
      subject: 'compound', code: 'EXC-1001', title: 'Clean hERG', finding: 'hERG IC50 22 uM, 1000x over potency'});
    const id = /\(([^)]+)\)$/.exec(report)![1];
    made.push(id);
    expect(report.startsWith('EXS-KIN1: observation recorded on compound EXC-1001 · Clean hERG'), true, report);
    const r = await retractObservation(id, 'the hERG run was repeated');
    expect(`${r.status} ${r.retraction}`, 'retracted the hERG run was repeated');
    await expectExceptionAsync(async () => {
      await retractObservation(id, 'again');
    });
    const feed = await programFeedOf(smd, kin1, {kinds: ['observation']});
    expect(feed.some((x) => x.text === 'Observation retracted on compound EXC-1001: Clean hERG — the hERG run was repeated'), true);
  });

  test('the dialog offers the program, its series and compounds', async () => {
    const dlg = await openObservationDialog({portfolio: smd, program: kin1});
    try {
      const options = Array.from(dlg.root.querySelectorAll('[data-u2-name="observation-subject-input"] option')).map((o) => o.textContent);
      expect(options[0], 'The program');
      expect(options.includes('series KIN-A') && options.includes('compound EXC-1001'), true, options.join());
      expect(!!dlg.root.querySelector('[data-u2-name="record-observation-button"]'), true);
    } finally {
      dlg.close();
    }
  });
});
