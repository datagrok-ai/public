import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {ADDED, applyPanels, quietPanels} from '../stages';
import type {Panel} from '../stages';
import {ncaInstalled} from '../exposure';
import {DEMO_TAG, demoBrought, demoStudies, hubView, openDemo, openStudyDialog, savedStudyPanel, showOpened}
  from '../hub-view';
import {helpPage, modal, modalOf, modalsHeld, preparing, showOnce, stageOf} from '../ui';
import {Dialog} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import {ONE_CLICK_MS, clickItem} from '../report';
import {readFile, routedNote, suggestedStage} from '../open-study';
import type {Brought} from '../open-study';
import {stageToolbox} from '../report';
import {analyzeView} from '../analyze-view';
import {monitorView} from '../monitor-view';
import type {SavedStudy} from '../study-file';

const CSV = 'Mouse ID,Arm,Day,Volume (mm3)\nm1,Vehicle,0,100\nm1,Vehicle,12,900\nm2,Cpd-1,0,105\nm2,Cpd-1,12,300\n';

/** A view as the panels read it: of a PMX stage (or none: another app's), pinned or not, added by the app or not. */
function viewOf(stage: string | null, pinned = false, added = false): DG.ViewBase {
  const root = document.createElement('div');
  if (stage)
    root.dataset.pmxStage = stage;
  if (added)
    root.dataset[ADDED] = 'true';
  return {root, isPinned: pinned} as unknown as DG.ViewBase;
}

/** The dialogs held once the closed are forgotten: one made and dropped (never shown) sweeps them. */
function openHeld(): number {
  modal('test-sweep');
  modalOf('test-sweep');
  return modalsHeld();
}

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

  test('the panels: Browse hidden only when the user pinned the app, given back as they were', async () => {
    const w: {[p in Panel]: boolean} = {showContextPanel: true, showBrowse: true};
    const before = new Map<Panel, boolean>();
    const app = {pinned: false};
    const panels = () => [w.showContextPanel, w.showBrowse].join();
    // the app opened as a preview (one click in Browse): the context panel hidden, Browse left
    expect(applyPanels(w, viewOf('home'), before, app), false);
    expect(panels(), 'false,true');
    // a stage the app added (the platform adds it pinned): it follows the app, Browse left
    expect(applyPanels(w, viewOf('plan', true, true), before, app), false);
    expect(panels(), 'false,true');
    // the app pinned (or opened by a double click): Browse hidden there and on the stages it opens
    expect(applyPanels(w, viewOf('home', true), before, app), true);
    expect(panels(), 'false,false');
    expect(applyPanels(w, viewOf('analyze', true, true), before, app), true);
    expect(panels(), 'false,false');
    // unpinned: Browse back
    applyPanels(w, viewOf('home'), before, app);
    expect(panels(), 'false,true');
    // leaving the app: as the user had it
    applyPanels(w, viewOf(null), before, app);
    expect(panels(), 'true,true');
    expect(before.size, 0);
    // Browse closed by the user: pinned and unpinned, it stays closed; the view closed (none current), all back
    w.showBrowse = false;
    applyPanels(w, viewOf('monitor', true), before, app);
    applyPanels(w, null, before, app);
    expect(panels(), 'true,false');
  });

  test('a dialog closed is forgotten: purposes made per view do not pile up', async () => {
    const held = openHeld();
    for (let i = 0; i < 3; i++) {
      const d = modal(`test-held-${i}`)!;
      d.show();
      d.close();
    }
    expect(openHeld(), held);
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

  test('a u2 dialog (Save plan as in Plan) clicked three times is one dialog; closed, it opens again', async () => {
    const view = grok.shell.addView(DG.View.create());
    const host = ViewHost.of(view);
    let made = 0;
    let shown: Dialog | null = null;
    const make = () => {
      made++;
      return shown = new Dialog('Save plan');
    };
    try {
      showOnce(host, 'test-u2', make);
      showOnce(host, 'test-u2', make);
      showOnce(host, 'test-u2', make);
      expect(made, 1);
      expect(shown!.isOpen.peek(), true);
      shown!.close('cancel');
      await delay(50);
      showOnce(host, 'test-u2', make);
      expect(made, 2);
    } finally {
      view.close();
    }
  });

  test('an export clicked three times fast runs once (one download)', async () => {
    let runs = 0;
    const item: [string, () => void] = ['Test export', () => void runs++];
    clickItem(item);
    clickItem(item);
    clickItem(item);
    await delay(50);
    expect(runs, 1);
    await delay(ONE_CLICK_MS + 50);
    clickItem(item);
    await delay(50);
    expect(runs, 2);
  });

  test('a dialog prepared (read before it shows): a second click meanwhile starts nothing', async () => {
    let runs = 0;
    let release = () => {};
    const slow = () => new Promise<void>((r) => {
      runs++;
      release = r;
    });
    const one = preparing('test-prepared', slow);
    expect(await preparing('test-prepared', slow) === undefined, true);
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
    // Plan's own Open…: still a study (Analyze), its use as a model's growth history offered and said
    const asHistory = await openedWithCsv({stage: 'plan', history: () => {}});
    expect(asHistory.ok, 'Open in Analyze');
    expect(asHistory.text.includes('growth history of a model outside the atlas'), true);
    const fromMonitor = await openedWithCsv({stage: 'monitor'});
    expect(fromMonitor.ok, 'Open in Analyze');
    expect(fromMonitor.text.includes('Without its plan'), true);
  });

  test('Open a study: one dialog at a time, a second click answered by the one open', async () => {
    const first = openStudyDialog();
    try {
      expect(openStudyDialog() === first, true);
      expect(document.querySelectorAll('.pmx-open').length, 1);
      // opened from a stage while the Overview's is open: the stage's takes its place
      const fromPlan = openStudyDialog(undefined, {stage: 'plan'});
      try {
        expect(fromPlan === first, false);
        expect(document.querySelectorAll('.pmx-open').length, 1);
        expect(openStudyDialog(undefined, {stage: 'plan'}) === fromPlan, true);
      } finally {
        fromPlan.close();
      }
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
    const monitors = () => [...grok.shell.views].filter((v) => stageOf(v) === 'monitor');
    const had = new Set(monitors().map((v) => v.root));
    await openDemo('monitor');
    // (the view may be added a moment after the demo is read)
    let opened = monitors().filter((v) => !had.has(v.root));
    for (let i = 0; i < 50 && !opened.length; i++) {
      await delay(100);
      opened = monitors().filter((v) => !had.has(v.root));
    }
    try {
      expect(opened.length, 1);
    } finally {
      opened.forEach((v) => v.close());
    }
  });

  test('the Analyze demo: its mouse PK opened in the workspace, for its exposure (with NCA)', async () => {
    if (!ncaInstalled())
      return;
    const analyzes = () => [...grok.shell.views].filter((v) => stageOf(v) === 'analyze');
    const had = new Set(analyzes().map((v) => v.root));
    await openDemo('analyze');
    let opened = analyzes().filter((v) => !had.has(v.root));
    for (let i = 0; i < 50 && !opened.length; i++) {
      await delay(100);
      opened = analyzes().filter((v) => !had.has(v.root));
    }
    try {
      expect(opened.length, 1);
      expect(grok.shell.tables.some((t) => t.name === 'ct26-trametinib-pk-mouse'), true);
    } finally {
      opened.forEach((v) => v.close());
      for (const t of grok.shell.tables.filter((x) => x.name.startsWith('ct26-trametinib')))
        grok.shell.closeTable(t);
    }
  });

  test('Compare two demo studies: saved once with their effect, found again the next time', async () => {
    const saved: SavedStudy[] = [];
    const first = await demoStudies([], async (doc) => void saved.push(doc));
    expect(saved.length, 2);
    expect(first.every((e) => e.kind === 'study' && e.stage === 3 && e.tags.includes(DEMO_TAG)), true);
    // each with the analysis's effect, so Compare reads the TGI rather than the median T/C
    expect(saved.every((d) => Object.keys(d.exposure?.effect?.tgi ?? {}).length > 0), true);
    // the same analysis Analyze runs: its evaluation day, its TGIs (rounded as shown), its control's rate and range
    const shown = (d: SavedStudy) => Object.values(d.exposure!.effect!.tgi).map((t) => Math.round(t))
      .sort((a, b) => a - b).join();
    const [ct26, mc38] = saved;
    expect(ct26.exposure!.effect!.day, 16);
    expect(shown(ct26), '12,24,63,87,90');
    expect(shown(mc38), '54,84');
    expect(saved.every((d) => Number.isFinite(d.exposure!.effect!.control?.slope ?? NaN)), true);
    expect(saved.every((d) => !!d.exposure!.effect!.ctd), true);
    // its panel (picked in Your studies): its best group from what the study kept, its facts, no owner's login
    const panel = savedStudyPanel(first[0], async () => JSON.stringify(ct26));
    for (let i = 0; i < 50 && !panel.querySelector('.pmx-object-verdict'); i++)
      await delay(20);
    const verdict = panel.querySelector('.pmx-object-verdict')?.textContent ?? '';
    expect(verdict.startsWith('Best group: '), true, verdict);
    expect(verdict.includes('TGI 90%'), true, verdict);
    expect(/[0-9] of 5 treated groups slow growth/.test(verdict), true, verdict);
    const text = panel.textContent ?? '';
    expect(['Model', 'Groups', 'Progress', 'Updated'].every((s) => text.includes(s)), true, text);
    // its groups as Analyze counts them (the control once); a finished study never "day 21 of 14"
    expect(text.includes(`Groups${ct26.exposure!.effect!.ctd!.groups.length},`), true, text);
    expect(ct26.exposure!.effect!.ctd!.groups.length, 6);
    expect(/Progress *Day [0-9]+, finished/.test(text), true, text);
    expect(text.includes(grok.shell.user.login), false, 'the owner’s login');
    // the second time: found among the user's, nothing saved again
    const again = await demoStudies(first, async (doc) => void saved.push(doc));
    expect(saved.length, 2);
    expect(again.map((e) => e.name).join(), first.map((e) => e.name).join());
    // one of them deleted (so not listed): only that one saved again
    const third = await demoStudies(first.slice(1), async (doc) => void saved.push(doc));
    expect(saved.length, 3);
    expect(third.length, 2);
  });

  test('the help panel of the app\'s views: the README published with the package', async () => {
    const view = await hubView();
    expect(view.helpUrl, helpPage());
    const r = await fetch(helpPage());
    expect(r.ok && (await r.text()).includes('## Opening a study'), true);
  });

  test('a view opened from Your studies: a second click brings it back, not a second view', async () => {
    const view = DG.View.create();
    view.root.dataset.pmxOpenedAs = 'test-opened';
    grok.shell.addView(view);
    try {
      expect(showOpened('test-opened'), true);
      expect(grok.shell.v.root === view.root, true);
      expect(showOpened('test-not-opened'), false);
    } finally {
      view.close();
    }
  });

  test('the demo files carry no one\'s login: owner, author and notes are "demo"', async () => {
    const plan = await demoBrought('plan');
    expect(plan.plan!.saved.author, 'demo');
    const monitor = await demoBrought('monitor');
    expect(monitor.study!.saved.meta?.owner, 'demo');
    expect(monitor.study!.saved.plan?.author ?? 'demo', 'demo');
  });

  test('the Overview toolbox: the ways into the app, the demo data and help', async () => {
    const view = await hubView();
    const text = view.toolbox?.textContent ?? '';
    for (const s of ['Open a study…', 'New plan', 'Compare saved studies', 'Compare two demo studies',
      'PMX Studio help'])
      expect(text.includes(s), true);
  });

  test('a stage’s toolbox: its menu’s items by section, the study’s links first, help last', async () => {
    let ran = '';
    const box = stageToolbox([['Report (PDF)…', () => {
      ran = 'report';
    }], ['ARRIVE checklist…', () => {}, 'Documents'], ['Layout…', () => {}, 'View'],
    ['Animals', () => {}, 'Data > Open as table'], ['Statistics', () => {}, 'Export > CSV']],
    [['Open…', () => {}]], [['PMX Studio help', () => {}]]);
    const panes = [...box.querySelectorAll('.d4-accordion-pane-header')].map((h) => h.textContent?.trim());
    expect(panes.join(), 'Study,Export,Documents,View,Help');
    const text = box.textContent ?? '';
    for (const s of ['Open…', 'CSV: Statistics', 'Open as table: Animals', 'ARRIVE checklist…'])
      expect(text.includes(s), true, s);
    [...box.querySelectorAll<HTMLElement>('a, .ui-link')].find((a) => a.textContent === 'Report (PDF)…')!.click();
    await delay(50);
    expect(ran, 'report');
  });

  test('Analyze has its toolbox: Open…, the report, its documents, the layout and help', async () => {
    const view = await analyzeView();
    try {
      const text = view.toolbox?.textContent ?? '';
      for (const s of ['Open…', 'Study details…', 'Report (PDF)…', 'ARRIVE checklist…', 'Layout…', 'PMX Studio help',
        'Study file (.pmxstudy)', 'Calendar (.ics)'])
        expect(text.includes(s), true, s);
    } finally {
      view.close();
    }
  });

  test('Monitor: forming the groups first, in the toolbox’s Study section', async () => {
    const view = await monitorView();
    try {
      // the Study section (after the platform's own panes): its items before the next section's
      const all = view.toolbox?.textContent ?? '';
      const text = all.slice(all.indexOf('Study'));
      expect(all.includes('Study'), true, all.slice(0, 80));
      for (const s of ['Open…', 'Propose groups…', 'Print randomization sheet'])
        expect(text.indexOf(s) >= 0 && text.indexOf(s) < text.indexOf('Export'), true, s);
    } finally {
      view.close();
    }
  });
});
