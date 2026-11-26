import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {quietPanels} from '../stages';
import * as grok from 'datagrok-api/grok';
import {demoBrought, hubView, openDemo, openStudyDialog} from '../hub-view';
import {modal, preparing} from '../ui';
import {readFile, routedNote, suggestedStage} from '../open-study';
import {stageOf} from '../ui';
import type {Brought} from '../open-study';

const CSV = 'Mouse ID,Arm,Day,Volume (mm3)\nm1,Vehicle,0,100\nm1,Vehicle,12,900\nm2,Cpd-1,0,105\nm2,Cpd-1,12,300\n';

/** The Open a study dialog opened from `at` with a CSV file, once it is read: its OK button and what it says. */
async function openedWithCsv(at: Parameters<typeof openStudyDialog>[1]): Promise<{ok: string, text: string}> {
  const d = openStudyDialog([new File([CSV], 'study.csv')], at);
  try {
    // (read with a short pause, so the change of state shows: see hub-view take)
    for (let i = 0; i < 40 && !d.root.querySelector('.pmx-read'); i++)
      await delay(100);
    return {ok: d.root.querySelector('.pmx-primary')?.textContent ?? '', text: d.root.textContent ?? ''};
  } finally {
    d.close();
  }
}

category('Hub: opening', () => {
  test('Browse is left as the user had it, and hidden only on a pinned view of the app', async () => {
    expect(quietPanels('home', false, false).join(), 'showContextPanel');
    expect(quietPanels('plan', false, false).join(), 'showContextPanel');
    expect(quietPanels(null, true, false).join(), 'showContextPanel');
    expect(quietPanels('home', false, true).join(), 'showContextPanel,showBrowse');
    expect(quietPanels(null, true, true).join(), 'showContextPanel,showBrowse');
    // not the app's view: nothing hidden, pinned or not
    expect(quietPanels(null, false, true).length, 0);
  });

  test('a dialog opened twice is one dialog; closed, it opens again', async () => {
    const first = modal('test-twice', 'Twice');
    first!.show();
    try {
      expect(modal('test-twice', 'Twice'), null);
    } finally {
      first!.close();
    }
    const again = modal('test-twice', 'Twice');
    expect(again !== null && again !== first, true);
    again!.close();
  });

  test('a dialog prepared (read before it shows): a second click meanwhile starts nothing', async () => {
    let runs = 0;
    let release = () => {};
    const slow = () => new Promise<void>((r) => {
      runs++;
      release = r;
    });
    const one = preparing('test-prepared', slow);
    expect(await preparing('test-prepared', slow), undefined);
    release();
    await one;
    expect(runs, 1);
    // done, it can be launched again
    const two = preparing('test-prepared', slow);
    release();
    await two;
    expect(runs, 2);
  });

  test('a file of another stage: where it opens, said, never refused', async () => {
    const table = readFile('study.csv', CSV);
    const b: Brought = table.kind === 'table' ? {table} : {};
    expect(routedNote(b, 'analyze', 'analyze'), null);
    expect(routedNote(b, 'plan', 'analyze'), 'A table of measurements opens in Analyze, not in Plan.');
    expect(routedNote(b, 'monitor', 'analyze')!.startsWith('Without its plan'), true);
    expect(routedNote({}, 'analyze', 'plan'), 'A plan opens in Plan, not in Analyze.');
  });

  test('Open a study from a stage: the table opens where it belongs, or in Plan as a growth history', async () => {
    const fromPlan = await openedWithCsv({stage: 'plan'});
    expect(fromPlan.ok, 'Open in Analyze');
    expect(fromPlan.text.includes('opens in Analyze, not in Plan'), true);
    // Plan's own Open…: the table may be a model's growth history, and that is what it opens as
    const asHistory = await openedWithCsv({stage: 'plan', history: () => {}});
    expect(asHistory.ok, 'Open in Plan');
    const fromMonitor = await openedWithCsv({stage: 'monitor'});
    expect(fromMonitor.ok, 'Open in Analyze');
    expect(fromMonitor.text.includes('Without its plan'), true);
  });

  test('Open a study: one dialog at a time, a second click answered by the one open', async () => {
    const first = openStudyDialog();
    try {
      expect(openStudyDialog() === first, true);
      expect(document.querySelectorAll('.pmx-open').length, 1);
    } finally {
      first.close();
    }
    // closed, the next click opens a new one
    const next = openStudyDialog();
    try {
      expect(next === first, false);
    } finally {
      next.close();
    }
  });

  test('the demo data: the study of each stage, read from the package, opens in that stage', async () => {
    const plan = await demoBrought('plan');
    expect(plan.plan?.name.length! > 0 && !plan.table && !plan.study, true);
    expect(suggestedStage(plan), 'plan');
    const monitor = await demoBrought('monitor');
    expect(!!monitor.study, true);
    expect(suggestedStage(monitor), 'monitor');
    // Analyze: the finished table with its plan, measured to its end
    const analyze = await demoBrought('analyze');
    expect(!!analyze.table && !!analyze.plan && !analyze.table.read.error, true);
    expect(analyze.table!.groups, 6);
    expect(suggestedStage(analyze), 'analyze');
  });

  test('the Monitor demo opens Monitor with the running study', async () => {
    const before = new Set([...grok.shell.views].map((v) => v.dart));
    await openDemo('monitor');
    const opened = [...grok.shell.views].filter((v) => !before.has(v.dart));
    try {
      expect(opened.some((v) => stageOf(v) === 'monitor'), true);
    } finally {
      opened.forEach((v) => v.close());
    }
  });

  test('the Overview toolbox: the ways into the app, the demo data and help', async () => {
    const view = await hubView();
    const text = view.toolbox?.textContent ?? '';
    for (const s of ['Open a study…', 'New plan', 'Compare saved studies', 'Compare two demo studies',
      'PMX Studio help'])
      expect(text.includes(s), true);
  });
});
