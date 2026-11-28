import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {ADDED, applyPanels, quietPanels} from '../stages';
import type {Panel} from '../stages';
import {ncaInstalled} from '../exposure';
import {DEMO_TAG, demoBrought, demoStudies, hubView, openBrought, openDemo, openStudyDialog, savedStudyPanel,
  showOpened}
  from '../hub-view';
import {helpPage, modal, modalOf, modalsHeld, preparing, showOnce, stageOf} from '../ui';
import {Dialog} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import {ONE_CLICK_MS, clickItem} from '../report';
import {readFile, routedNote, suggestedStage} from '../open-study';
import type {Brought, PlanRead} from '../open-study';
import {_package} from '../package';
import {fileEntries, reasonOf, relabel, stageToolbox} from '../report';
import {crumb, headBar, ribbonEnd, steps} from '../study-header';
import {pageIndex, pageOf, panelGrid} from '../panels';
import {fitSpan, insertionIndex, snapSpan, stepSpan} from '../layout-motion';
import {openStage} from '../stages';
import {analyzeView} from '../analyze-view';
import {monitorView} from '../monitor-view';
import {newStudyTable} from '../study-file';
import type {SavedStudy} from '../study-file';
import {parsePlan} from '../plan-file';
import {emptyMeta} from '../study-meta';
import {appTable} from '../native';

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
    expect(quietPanels('plan', false, false).join(), '', 'a study\'s stage keeps the context panel: its study');
    expect(quietPanels('analyze', false, true).join(), 'showBrowse');
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
    // a stage the app added (the platform adds it pinned): it follows the app, Browse left, the context panel back
    expect(applyPanels(w, viewOf('plan', true, true), before, app), false);
    expect(panels(), 'true,true');
    // the app pinned (or opened by a double click): Browse hidden there and on the stages it opens
    expect(applyPanels(w, viewOf('home', true), before, app), true);
    expect(panels(), 'false,false');
    expect(applyPanels(w, viewOf('analyze', true, true), before, app), true);
    expect(panels(), 'true,false');
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

  test('a plan under a .csv name opens as a plan; a converted volume is not also "as given" (QA round 4)', async () => {
    const plan = await _package.files.readAsText('demo/ct26-trametinib-antipd1.pmxplan');
    const r = readFile('renamed.csv', plan);
    expect([r.kind, r.file].join(), 'plan,renamed.csv');
    const cm = readFile('cm.csv', 'animal,group,day,volume (cm3)\nm1,A,0,0.1\nm1,A,3,0.2\nm2,B,0,0.1\nm2,B,3,0.3\n');
    const notes = cm.kind === 'table' ? cm.read.notes : [];
    expect([notes.includes('Volumes as given.'), notes.some((n) => n.includes('converted'))].join(), 'false,true',
      notes.join(' | '));
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
      // the Overview's offers the saved studies as Monitor's Ctrl+O does; Plan's only the plans (QA 0.1.3 r5)
      const links = (d: DG.Dialog) => [...d.root.querySelectorAll('.pmx-alts a, .pmx-alts .ui-link')]
        .map((a) => a.textContent?.trim()).join();
      expect(links(first), 'Saved plans,Saved studies');
      // opened from a stage while the Overview's is open: the stage's takes its place
      const fromPlan = openStudyDialog(undefined, {stage: 'plan'});
      try {
        expect(links(fromPlan), 'Saved plans');
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
      // the app's own table: the untouched demo leaves no platform "Leave site?" (QA 0.1.3 r4)
      expect(appTable(grok.shell.tables.find((t) => t.name === 'ct26-trametinib-pk-mouse')!), true);
    } finally {
      opened.forEach((v) => v.close());
      for (const t of grok.shell.tables.filter((x) => x.name.startsWith('ct26-trametinib')))
        grok.shell.closeTable(t);
    }
  });

  test('a plan opened in Monitor whose study is followed opens that study, not a second unsaved one (QA 0.1.3 r4)',
    async () => {
      const plan = parsePlan(await _package.files.readAsText('demo/ct26-trametinib-antipd1.pmxplan'));
      const withId = {...plan, study: {...plan.study ?? emptyMeta(), id: 'PMX-TEST-FOLLOWED'}};
      const followed = newStudyTable(withId);
      grok.shell.addTable(followed);
      const monitors = () => [...grok.shell.views].filter((v) => stageOf(v) === 'monitor');
      const had = new Set(monitors().map((v) => v.root));
      const tables = grok.shell.tables.length;
      try {
        await openBrought({plan: readFile('followed.pmxplan', JSON.stringify(withId)) as PlanRead}, 'monitor');
        await delay(500);
        expect(grok.shell.tables.length, tables, 'no second table of the same study');
        expect(monitors().filter((v) => !had.has(v.root)).length, 1, 'Monitor on the followed study');
      } finally {
        monitors().filter((v) => !had.has(v.root)).forEach((v) => v.close());
        grok.shell.closeTable(followed);
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
    expect(r.ok && (await r.text()).includes('## Try it'), true);
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

  test('a stage’s toolbox: the page’s index, then the actions not in File, help last', async () => {
    let ran = '';
    const box = stageToolbox([['Report (PDF)…', () => {}], ['ARRIVE checklist…', () => {}, 'Documents'],
      ['Customize this page', () => {
        ran = 'customize';
      }, 'View'], ['Propose groups…', () => {}, 'Study'], ['Statistics', () => {}, 'Export > CSV']],
    [['Open…', () => {}], ['Study details…', () => {}]], [['PMX Studio help', () => {}]]);
    const panes = [...box.querySelectorAll('.d4-accordion-pane-header')].map((h) => h.textContent?.trim());
    expect(panes.join(), 'Actions');
    const text = box.textContent ?? '';
    for (const s of ['Study details…', 'Propose groups…', 'Customize this page', 'PMX Studio help'])
      expect(text.includes(s), true, s);
    // Open, Save and the exports are File's
    for (const s of ['Open…', 'Report (PDF)…', 'ARRIVE checklist…', 'Statistics', 'in File'])
      expect(text.includes(s), false, s);
    [...box.querySelectorAll<HTMLElement>('a, .ui-link')].find((a) => a.textContent === 'Customize this page')!.click();
    await delay(50);
    expect(ran, 'customize');
  });

  test('File: Open and Save with their keys, then Export, Documents and View; an item off says why', async () => {
    const entries = fileEntries([['Report (PDF)…', () => {}], ['Save as…', () => {}, 'File', () => 'Plan first'],
      ['ARRIVE checklist…', () => {}, 'Documents'], ['Customize this page', () => {}, 'View'],
      ['Animals', () => {}, 'Data > Open as table'], ['Statistics', () => {}, 'Export > CSV']],
    {open: () => {}, save: {text: () => 'Save plan', run: () => {}, why: () => 'Add a treatment arm first'}});
    // (separators between the groups)
    type Item = {text: string, why: () => string | null, shortcut?: string};
    const items = entries.filter((e) => e !== '-') as Item[];
    expect(items.map((e) => e.text).join(), 'Open…,Save plan,Save as…,Export,Documents,View');
    const [open, save, saveAs] = items;
    expect(open.shortcut, 'Ctrl+O');
    expect(open.why(), null);
    expect(save.shortcut, 'Ctrl+S');
    expect(save.why(), 'Add a treatment arm first');
    expect(saveAs.why(), 'Plan first');
    const view = entries.find((e) => e !== '-' && e.text === 'View') as {sub: {text: string}[]};
    expect(JSON.stringify(view.sub).includes('Customize this page') &&
      JSON.stringify(view.sub).includes('Animals'), true);
    // an item of a group File and the toolbox do not know is File's, not lost; the toolbox's Study items are not
    const odd = fileEntries([['Odd one', () => {}, 'Elsewhere'], ['Propose groups…', () => {}, 'Study']],
      {open: () => {}}).filter((e) => e !== '-') as Item[];
    expect(odd.map((e) => e.text).join(), 'Open…,Odd one');
    // a button hidden or off gives its reason
    const b = document.createElement('button');
    expect(reasonOf(b), null);
    b.hidden = true;
    expect(reasonOf(b), 'Open a study first');
  });

  // QA 0.1.3 browser check: the ribbon said "Open the followed study", the toolbox "Follow this study"
  test('an item named as its button: the toolbox and File follow the name as it changes', async () => {
    let name = 'Follow this study';
    const items: Parameters<typeof stageToolbox>[0] = [[() => name, () => {}, 'Study'],
      [() => name, () => {}, 'Export']];
    const box = stageToolbox(items, [], []);
    const text = () => box.textContent ?? '';
    expect(text().includes('Follow this study'), true);
    name = 'Open the followed study';
    relabel(box);
    expect([text().includes('Open the followed study'), text().includes('Follow this study')].join(), 'true,false');
    expect(JSON.stringify(fileEntries(items, {open: () => {}})).includes('Open the followed study'), true);
  });

  // QA 0.1.3 browser check: Tab went Discard, Save, Cancel over Cancel · Save · Discard
  test('every dialog: its buttons in the page in the order they show', async () => {
    const d = modal('test-button-order', 'Unsaved changes')!;
    d.add(document.createElement('div')).onOK(() => {});
    d.addButton('Discard', () => d.close());
    d.show();
    try {
      const buttons = [...d.root.querySelectorAll<HTMLElement>('.d4-dialog-footer button, .d4-command-bar button')]
        .filter((b) => b.offsetParent !== null);
      expect(buttons.length >= 3, true, 'OK, Cancel and Discard');
      const lefts = buttons.map((b) => b.getBoundingClientRect().left);
      expect(lefts.every((x, i) => i === 0 || x > lefts[i - 1]), true,
        `page order ${buttons.map((b) => b.textContent?.trim()).join(', ')}`);
    } finally {
      d.close();
    }
  });

  test('the home icon before the steps leads to the Overview; the steps are the three stages, no Compare', async () => {
    let to = '';
    const c = crumb(() => to = 'home');
    // (no text: no "PMX Studio › Plan", the steps say the stage)
    expect(c.textContent, '');
    expect(c.querySelector('.pmx-crumb-home')?.getAttribute('aria-label'), 'Overview');
    c.querySelector<HTMLElement>('.pmx-crumb-home')!.click();
    expect(to, 'home');
    expect(crumb(() => {}, 'Compare · 2 studies').textContent, 'Compare · 2 studies');
    const nav = steps('analyze', () => {});
    expect([...nav.querySelectorAll<HTMLElement>('.pmx-head-step')].map((s) => s.dataset.stage).join(),
      'plan,monitor,analyze');
  });

  // QA 0.1.3 browser check: Plan at 1600 px, the ⌂ and Plan under the ribbon's Study and Help (elementFromPoint(560,
  // 16)): the menu's box ended before its items did
  test('the stage bar starts after the ribbon\'s items, measured, and again as they change', async () => {
    const box = (cls: string, style: Partial<CSSStyleDeclaration> = {}) => {
      const e = document.createElement('div');
      e.className = cls;
      Object.assign(e.style, style);
      return e;
    };
    const row = box('d4-ribbon', {position: 'fixed', left: '0', top: '0', width: '1600px', height: '30px',
      display: 'flex', zIndex: '-1'});
    // the menu: a narrow box whose items run past it
    const menu = box('', {width: '60px', overflow: 'visible', whiteSpace: 'nowrap', display: 'flex', flex: 'none'});
    menu.setAttribute('role', 'menubar');
    const item = (t: string) => Object.assign(document.createElement('span'), {textContent: t,
      style: 'display:inline-block;width:120px'});
    menu.append(item('File'), item('View'));
    const inner = box('d4-ribbon', {flex: '1 1 0', display: 'flex', minWidth: '0'});
    const center = box('', {width: '200px'});
    center.append(item('⌂ Plan'));
    const bar = headBar(center, box(''));
    inner.append(bar);
    row.append(menu, inner);
    document.body.append(row);
    try {
      const frames = async () => {
        for (let i = 0; i < 4; i++)
          await new Promise((r) => requestAnimationFrame(r));
      };
      await frames();
      expect(ribbonEnd(row, menu, bar), menu.lastElementChild!.getBoundingClientRect().right);
      const clear = () => center.getBoundingClientRect().left >= ribbonEnd(row, menu, bar);
      expect(clear(), true, 'the steps after File and View');
      // more items put in, no size of the menu's box changed: placed again
      menu.append(item('Study'), item('Help'));
      await frames();
      expect(clear(), true, 'the steps after Study and Help');
    } finally {
      row.remove();
    }
  });

  test('a folded panel stays folded for the user; Reset gives the stage\'s folds back', async () => {
    const stage = 'test-fold';
    grok.userSettings.delete('PmxStudio', `fold-${stage}`);
    grok.userSettings.delete('PmxStudio', `layout-${stage}`);
    const panel = (id: string) => ({id, title: id, hint: '', content: [document.createElement('div')]});
    const layout = [{id: 'a', span: 12 as const}, {id: 'b', span: 12 as const, folded: true}];
    try {
      const page = panelGrid([panel('a'), panel('b')], layout, 'test', stage);
      const cell = (root: HTMLElement, id: string) => root.querySelector<HTMLElement>(`[data-panel="${id}"]`)!;
      expect(cell(page, 'b').classList.contains('is-folded'), true);
      expect(cell(page, 'a').classList.contains('pmx-foldable'), false);
      cell(page, 'b').querySelector<HTMLElement>('.pmx-block-header')!.click();
      expect(cell(page, 'b').classList.contains('is-folded'), false);
      // another page of the stage (the next study) opens it unfolded
      const again = panelGrid([panel('a'), panel('b')], layout, 'test', stage);
      expect(cell(again, 'b').classList.contains('is-folded'), false);
      pageOf(again)!.reset();
      expect(cell(again, 'b').classList.contains('is-folded'), true);
      expect(panelGrid([panel('a'), panel('b')], layout, 'test', stage).querySelector('[data-panel="b"]')!
        .classList.contains('is-folded'), true);
    } finally {
      grok.userSettings.delete('PmxStudio', `fold-${stage}`);
      grok.userSettings.delete('PmxStudio', `layout-${stage}`);
    }
  });

  test('customizing: a drag lands in the nearest slot, an edge snaps to a width, Ctrl+Z takes a change back',
    async () => {
      // two panels side by side on a row (0-400, 420-820), a third under them (0-820, 300 down)
      const boxes = [{left: 0, top: 0, width: 400, height: 200}, {left: 420, top: 0, width: 400, height: 200},
        {left: 0, top: 300, width: 820, height: 200}];
      expect(insertionIndex(boxes, 100, 50), 0, 'over the first, its left half: before it');
      expect(insertionIndex(boxes, 300, 50), 1, 'its right half: after it');
      expect(insertionIndex(boxes, 410, 100), 1, 'in the gap between the two: between them');
      expect(insertionIndex(boxes, 700, 260), 2, 'between the rows, nearer the third: before it');
      expect(insertionIndex(boxes, 400, 900), 3, 'far below everything: last');
      expect(insertionIndex(boxes, -500, -500), 0, 'far above and left: first');
      expect(insertionIndex([], 10, 10), 0);
      expect([2, 3.4, 3.6, 5.1, 6.9, 7.1, 9.9, 10.1, 20].map(snapSpan).join(), '3,3,4,6,6,8,8,12,12');
      // ← / → on an edge: one stop, held at a quarter and the row
      expect([stepSpan(3, -1), stepSpan(3, 1), stepSpan(6, 1), stepSpan(8, -1), stepSpan(12, 1), stepSpan(5, 1)].join(),
        '3,4,8,6,12,6');
      // dragged into a row with a gap, a panel takes the gap; where it fits, or there is no gap, its own width
      expect(fitSpan([6], 12), 6, 'a whole row beside a half: a half');
      expect(fitSpan([4], 12), 8, 'beside a third: two thirds');
      expect(fitSpan([3, 3], 8), 6, 'beside two quarters: a half');
      expect(fitSpan([6], 4), 4, 'fits: its own');
      expect(fitSpan([12], 6), 6, 'a full row, no gap: its own (a row of its own)');
      expect(fitSpan([6, 4], 12), 12, 'a gap narrower than a quarter: its own');
      expect(fitSpan([], 12), 12);
      const stage = 'test-undo';
      grok.userSettings.delete('PmxStudio', `layout-${stage}`);
      const panel = (id: string) => ({id, title: id, hint: '', content: [document.createElement('div')]});
      const root = panelGrid([panel('a'), panel('b')], [{id: 'a', span: 12}, {id: 'b', span: 12}], 'test', stage);
      document.body.append(root);
      try {
        const page = pageOf(root)!;
        const cell = (id: string) => root.querySelector<HTMLElement>(`[data-panel="${id}"]`);
        const shown = () => [...root.querySelectorAll<HTMLElement>('.pmx-panel')].map((c) =>
          `${c.dataset.panel}:${c.dataset.span}`).join();
        page.customize(true);
        const edge = cell('a')!.querySelector<HTMLElement>('.pmx-resize')!;
        edge.dispatchEvent(new KeyboardEvent('keydown', {key: 'ArrowLeft', bubbles: true}));
        edge.dispatchEvent(new KeyboardEvent('keydown', {key: 'ArrowLeft', bubbles: true}));
        expect(shown(), 'a:6,b:12');
        expect(edge.getAttribute('aria-valuetext'), 'Half the row');
        page.undo();
        expect(shown(), 'a:8,b:12', 'each step its own undo');
        edge.dispatchEvent(new KeyboardEvent('keydown', {key: 'ArrowLeft', bubbles: true}));
        cell('b')!.querySelector<HTMLElement>('.pmx-hide')!.click();
        await delay(400); // (it fades out first)
        expect(shown(), 'a:6');
        cell('a')!.querySelector<HTMLElement>('.pmx-handle')!.dispatchEvent(new KeyboardEvent('keydown',
          {key: 'ArrowDown', bubbles: true}));
        // Ctrl+Z, newest first: the move (none: one panel shown, it stays), the hide, the width
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'z', ctrlKey: true, bubbles: true}));
        expect(shown(), 'a:6');
        expect(page.undo(), true);
        expect(shown(), 'a:6,b:12');
        expect(page.undo(), true);
        expect(shown(), 'a:8,b:12');
        expect(page.undo(), true);
        expect(shown(), 'a:12,b:12');
        expect(page.undo(), false, 'nothing left to take back');
        page.cancel();
      } finally {
        root.remove();
        grok.userSettings.delete('PmxStudio', `layout-${stage}`);
      }
    });

  // QA 0.1.3 browser check: Customize with the context panel open, its instructions one word a line
  test('Customize on a narrow page: its instructions wrap as text, a line of their own if need be', async () => {
    const stage = 'test-narrow-bar';
    grok.userSettings.delete('PmxStudio', `layout-${stage}`);
    const panel = (id: string) => ({id, title: id, hint: '', content: [document.createElement('div')]});
    const root = panelGrid([panel('a')], [{id: 'a', span: 12}], 'test', stage);
    const host = document.createElement('div');
    host.style.width = '480px';
    host.append(root);
    document.body.append(host);
    try {
      pageOf(root)!.customize(true);
      const help = root.querySelector<HTMLElement>('.pmx-bar-text')!;
      expect(help.getBoundingClientRect().width >= 200, true, `${help.getBoundingClientRect().width} px wide`);
      pageOf(root)!.cancel();
    } finally {
      host.remove();
      grok.userSettings.delete('PmxStudio', `layout-${stage}`);
    }
  });

  test('the page’s index lists the hidden panels greyed, with + show', async () => {
    const stage = 'test-index';
    grok.userSettings.delete('PmxStudio', `layout-${stage}`);
    const panel = (id: string) => ({id, title: `Panel ${id}`, hint: '', content: [document.createElement('div')]});
    try {
      const page = panelGrid([panel('a'), panel('b')], [{id: 'a', span: 12}, {id: 'b', span: 12, hidden: true}],
        'test', stage);
      const index = pageIndex(page);
      const hidden = index.querySelector<HTMLElement>('.pmx-index-item.is-hidden')!;
      expect(hidden.textContent?.includes('Panel b') && hidden.textContent.includes('+ show'), true);
      expect((page.querySelector('.pmx-hidden-line') as HTMLElement).hidden, false);
      hidden.querySelector<HTMLElement>('.pmx-index-show')!.click();
      await delay(20);
      expect(index.querySelector('.pmx-index-item.is-hidden'), null);
      expect(page.querySelector('[data-panel="b"]') !== null, true);
    } finally {
      grok.userSettings.delete('PmxStudio', `layout-${stage}`);
    }
  });

  test('Analyze: the page’s index and its actions in the toolbox, File in the view menu', async () => {
    const view = await analyzeView();
    try {
      const text = view.toolbox?.textContent ?? '';
      for (const s of ['On this page', 'Experimental groups', 'Study details…', 'Customize this page',
        'PMX Studio help'])
        expect(text.includes(s), true, s);
      // the platform's own view menu, File one of its groups
      const menu = view.ribbonMenu.root.textContent ?? '';
      for (const g of ['File', 'View', 'Study', 'Help'])
        expect(menu.includes(g), true, g);
    } finally {
      view.close();
      ViewHost.of(view).dispose();
    }
  });

  test('Monitor: forming the groups first, among the toolbox’s actions', async () => {
    const view = await monitorView();
    try {
      const all = view.toolbox?.textContent ?? '';
      const text = all.slice(all.indexOf('Actions'));
      expect(all.includes('Actions'), true, all.slice(0, 80));
      for (const s of ['Propose groups…', 'Print randomization sheet'])
        expect(text.indexOf(s) >= 0 && text.indexOf(s) < text.indexOf('PMX Studio help'), true, s);
    } finally {
      view.close();
      ViewHost.of(view).dispose();
    }
  });

  test('one view: moving between the stages replaces the one shown', async () => {
    const stageViews = () => [...grok.shell.views].filter((v) => ['plan', 'monitor', 'analyze']
      .includes(stageOf(v) ?? ''));
    const before = stageViews().length;
    try {
      await openStage('analyze');
      await delay(300);
      expect(stageViews().filter((v) => stageOf(v) === 'analyze').length, 1);
      await openStage('monitor');
      await delay(500);
      expect(stageViews().filter((v) => stageOf(v) === 'monitor').length, 1);
      expect(stageViews().filter((v) => stageOf(v) === 'analyze').length, 0);
      expect(stageViews().length <= Math.max(1, before), true);
    } finally {
      stageViews().forEach((v) => v.close());
    }
  });
});
