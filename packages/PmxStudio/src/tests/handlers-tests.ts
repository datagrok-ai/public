import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {GROUP_ACTIONS, PmxAnimal, PmxFigure, PmxGroup, actionsOf, figurePanel, registerHandlers, withActions}
  from '../handlers';
import {animalMarks, mark} from '../study-table';
import {AnimalPanels, AnimalView, growthFrame, liveAnimalCharts, ordinal, placeWords, rankOf, weightFrame}
  from '../animal-panel';
import {userSetting, setUserSetting, sized} from '../ui';
import {EXCLUDED} from '../study-file';
import type {StudyInputs} from '../explore';

// Animals, groups and figures as the platform's objects: the handlers take them, each action offered when it
// applies, and the actions do what the old menus did (a mark cleared is the same mark cleared).
category('Right-click: animals, groups and figures', () => {
  const study = (): StudyInputs => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M2', 'M2']),
      DG.Column.fromStrings('group', ['Vehicle', 'Vehicle', 'Cpd-A', 'Cpd-A']),
      DG.Column.fromList(DG.TYPE.INT, 'day', [0, 7, 0, 7]),
      DG.Column.fromList(DG.TYPE.FLOAT, 'volume', [100, 300, 100, 150])]);
    return {study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null} as
      StudyInputs;
  };
  const names = (xs: {name: string}[]) => xs.map((a) => a.name);

  test('each handler takes its own objects, and only those', async () => {
    registerHandlers();
    const handler = (type: string) => DG.ObjectHandler.list().find((h) => h.type === type)!;
    const [animal, group, figure] = ['pmx-animal', 'pmx-group', 'pmx-figure'].map(handler);
    const n = study();
    expect(animal.isApplicable(new PmxAnimal(n, 'M1')), true);
    expect(animal.isApplicable(new PmxGroup(n, 'Vehicle', '#888')), false);
    expect(group.isApplicable(new PmxGroup(n, 'Vehicle', '#888')), true);
    expect(figure.isApplicable(new PmxFigure('growth')), true);
    expect(figure.isApplicable('growth'), false);
    // no tooltip of the platform's: the element's own is the one shown; the app's help in the help panel
    expect([animal, group, figure].every((h) => (h.helpUrl ?? '').endsWith('README.md')), true);
    expect([animal, group, figure].every((h) => h.renderTooltip(new PmxGroup(n, 'A', '#888')) === null), true);
    // the details shown on a click, then an Actions section with every action the right click offers
    const props = animal.renderProperties(new PmxAnimal(n, 'M1'));
    // (but those of the right-click menu only)
    const shown = (xs: {name: string, menuOnly?: boolean}[]) => names(xs.filter((a) => !a.menuOnly));
    expect(shown(actionsOf(new PmxAnimal(n, 'M1'))).every((a) => props.textContent?.includes(a)), true,
      props.textContent ?? '');
    expect(props.textContent?.includes('Actions'), true);
    const g = group.renderProperties(new PmxGroup(n, 'Vehicle', '#888'));
    expect(shown(actionsOf(new PmxGroup(n, 'Vehicle', '#888'))).every((a) => g.textContent?.includes(a)), true,
      g.textContent ?? '');
  });

  test('an animal: no "Show against its twin" (its panel draws it); Exclude, or Include once excluded', async () => {
    const n = study();
    const plain = names(actionsOf(new PmxAnimal(n, 'M1')));
    expect(plain.includes('Show against its twin'), false);
    expect(plain.includes('Exclude from the analysis…') && !plain.includes('Include in the analysis'), true);
    // excluded: Include is offered, and runs as the old Clear its marks did
    mark(n.study, [2, 3], 'excluded', 'tumor ulcerated', 'animal M2 (2 rows)');
    const after = names(actionsOf(new PmxAnimal(n, 'M2')));
    expect(after.includes('Include in the analysis') && !after.includes('Exclude from the analysis…'), true,
      after.join());
    actionsOf(new PmxAnimal(n, 'M2')).find((a) => a.name === 'Include in the analysis')!.run(new PmxAnimal(n, 'M2'));
    expect(animalMarks(n).get('M2') === undefined, true);
    expect(n.study.col(EXCLUDED)!.get(2), false);
  });

  // one animal panel, whichever chart or list the animal is clicked in (Analyze's treated and control, Monitor's)
  const weighed = (): StudyInputs => {
    const n = study();
    n.study.columns.add(DG.Column.fromList(DG.TYPE.FLOAT, 'weight', [20, 21, 20, 18.8]));
    return {...n, bodyWeight: 'weight'} as StudyInputs;
  };
  const treatedView = (animal: string): AnimalView => animal === 'M2' ? {
    verdict: {text: 'Responder: 40% smaller than its twin', tone: 'ok'},
    bar: Object.assign(document.createElement('div'), {className: 'pmx-animal-effect'}),
    curve: {title: 'Against its twin', frame: growthFrame([0, 7], [100, 150], '#d62728')!},
    group: {bars: [{animal: 'M2', group: 'Cpd-A', value: -40}, {animal: 'M3', group: 'Cpd-A', value: 10}],
      order: 'desc', summary: '1st of 2, −40% best'}, weightLimit: 20} :
    {curve: {title: 'Growth', frame: growthFrame([0, 7], [100, 300], '#1f77b4')!}, weightLimit: 20,
      group: {bars: [{animal: 'M1', group: 'Vehicle', value: 200}], order: 'desc', summary: '1st of 1, +200% best'}};
  const panelsOf = (n: StudyInputs, o: {color?: (a: string) => string | undefined,
    view?: (a: string) => AnimalView, onShow?: (a: string) => void} = {}) => new AnimalPanels({names: () => n,
    color: o.color ?? ((a) => a === 'M2' ? '#d62728' : '#1f77b4'), view: o.view ?? treatedView,
    onShow: o.onShow});
  // (the pane names in order, and which are open)
  const panesOf = (p: AnimalPanels) => ((p as any).panes as [string, DG.AccordionPane][])
    .map(([t, pane]) => `${t}${pane.expanded ? '+' : ''}`).join(', ');
  // (the user's kept folds put back after each test)
  const keepingFolds = async (body: () => Promise<void> | void) => {
    const was = userSetting('animalPanelFolds');
    setUserSetting('animalPanelFolds', '{}');
    try {
      await body();
    } finally {
      setUserSetting('animalPanelFolds', was ?? '{}');
    }
  };

  test('a treated animal: name and group, the answer, its bar under it; its twin open, the rest folded', async () =>
    keepingFolds(() => {
      const n = weighed();
      const p = panelsOf(n);
      try {
        const el = p.render('M2');
        expect(el.querySelector('.pmx-object-title')?.textContent, 'M2');
        expect(el.querySelector('.pmx-object-subtitle')?.textContent, 'Cpd-A');
        expect(el.querySelector('.pmx-object-verdict')?.textContent, 'Responder: 40% smaller than its twin');
        expect(el.querySelector('.pmx-object-verdict')?.nextElementSibling?.className, 'pmx-animal-effect');
        // no "90% range … to …" line: the bar says it
        expect(el.textContent?.includes('90% range'), false);
        expect(panesOf(p), 'Against its twin+, Among its group, Body weight, Details');
      } finally {
        p.dispose();
      }
    }));

  test('a control: no answer nor bar, its growth instead of a twin; Monitor\'s, blinded: no group', async () =>
    keepingFolds(() => {
      const n = weighed();
      const p = panelsOf(n);
      const blind = panelsOf(n, {color: () => undefined,
        view: () => ({curve: {title: 'Growth', frame: growthFrame([0, 7], [100, 300], undefined,
          {day: 10, lower: 300, median: 400, upper: 520})!}, weightLimit: 20})});
      try {
        const el = p.render('M1');
        expect(el.querySelector('.pmx-object-verdict'), null);
        expect(el.querySelector('.pmx-animal-effect'), null);
        expect(panesOf(p), 'Growth+, Among its group, Body weight, Details');
        const m = blind.render('M1');
        expect(m.querySelector('.pmx-object-subtitle'), null, 'a blinded study hides the group');
        expect(m.querySelector('.pmx-object-verdict'), null);
        expect(panesOf(blind), 'Growth+, Body weight, Details');
        // without body weight: no Body weight section
        const plain = panelsOf(study());
        plain.render('M1');
        expect(panesOf(plain), 'Growth+, Among its group, Details');
        plain.dispose();
      } finally {
        p.dispose();
        blind.dispose();
      }
    }));

  test('the same panel from two entry points (a row\'s object, and show)', async () => keepingFolds(async () => {
    const n = weighed();
    const shown: string[] = [];
    const p = panelsOf(n, {onShow: (a) => shown.push(a)});
    registerHandlers();
    const handler = DG.ObjectHandler.list().find((h) => h.type === 'pmx-animal')!;
    try {
      const fromRow = handler.renderProperties(p.objectOf('M2')).textContent;
      p.show('M2');
      await delay(50);
      const now = grok.shell.o as PmxAnimal;
      expect(now instanceof PmxAnimal && now.animal === 'M2', true);
      expect(handler.renderProperties(now).textContent, fromRow);
      expect(shown.join(), 'M2', 'marked on the page');
    } finally {
      grok.shell.o = null;
      p.dispose();
    }
  }));

  test('the folds are kept from one animal to the next', async () => keepingFolds(() => {
    const n = weighed();
    const p = panelsOf(n);
    try {
      p.render('M2');
      const panes = (p as any).panes as [string, DG.AccordionPane][];
      panes.find(([t]) => t === 'Details')![1].expanded = true;
      panes.find(([t]) => t === 'Against its twin')![1].expanded = false;
      p.render('M1');
      expect(panesOf(p), 'Growth+, Among its group, Body weight, Details+');
      p.render('M2');
      expect(panesOf(p), 'Against its twin, Among its group, Body weight, Details+');
      expect(JSON.parse(userSetting('animalPanelFolds') ?? '{}').Details, true, 'kept per user');
    } finally {
      p.dispose();
    }
  }));

  test('no chart left behind: five animals, three charts at most; none once something else is shown', async () =>
    keepingFolds(async () => {
      setUserSetting('animalPanelFolds', JSON.stringify({'Against its twin': true, 'Growth': true,
        'Among its group': true, 'Body weight': true}));
      const n = weighed();
      const p = panelsOf(n);
      const before = liveAnimalCharts();
      try {
        for (const a of ['M1', 'M2', 'M1', 'M2', 'M1'])
          document.body.append(p.render(a));
        expect(liveAnimalCharts() - before <= 3, true, String(liveAnimalCharts() - before));
        grok.shell.o = new PmxGroup(n, 'Vehicle', '#888');
        await delay(50);
        expect(liveAnimalCharts(), before, 'another object current lets them go');
      } finally {
        document.querySelectorAll('.pmx-object-panel').forEach((x) => x.remove());
        grok.shell.o = null;
        p.dispose();
      }
    }));

  test('an animal\'s place in its group: smallest first, in words', async () => {
    const bars = [{animal: 'a', group: 'G', value: -10}, {animal: 'b', group: 'G', value: -42},
      {animal: 'c', group: 'G', value: 5}];
    expect(rankOf(bars, 'b'), 1);
    expect(rankOf(bars, 'c'), 3);
    expect(rankOf(bars, 'x'), 0);
    expect([1, 2, 3, 4, 11, 12, 13, 21, 22, 101].map(ordinal).join(), '1st,2nd,3rd,4th,11th,12th,13th,21st,22nd,101st');
    // counted from the nearer end: the largest is "the largest", not "10th smallest of 10" (A11)
    const ten = Array.from({length: 10}, (_, i) => ({animal: `A${i + 1}`, group: 'G', value: i}));
    expect(['A1', 'A3', 'A9', 'A10', 'x'].map((a) => placeWords(ten, a)).join('|'),
      'the smallest of 10|3rd smallest of 10|2nd largest of 10|the largest of 10|');
    // its weight in words: the worst loss and its day; one that rounds to 0% is none
    expect(weightFrame([[0, 20], [2, 18.8], [4, 19.5]], undefined)!.summary, 'Lost at most 6.0% (day 2); day 4: −2%');
    // to the tenth as the tables, with the limit (designer round 3 #4: "−20% at most" beside a table's 19.6%)
    expect(weightFrame([[0, 20], [9, 16.08]], undefined, 20)!.summary, 'Lost at most 19.6% (day 9, limit 20%)');
    expect(weightFrame([[0, 20], [2, 19.95]], undefined)!.summary, 'No loss');
  });

  test('a panel chart takes the pane width in a flex pane (it was 0 px wide: A3, T10)', async () => {
    const pane = document.createElement('div');
    Object.assign(pane.style, {display: 'flex', width: '300px'});
    const root = document.createElement('div');
    root.style.position = 'relative';
    root.append(Object.assign(document.createElement('div'), {style: 'position: absolute; inset: 0'}));
    pane.append(sized({root}, 80));
    document.body.append(pane);
    try {
      expect(root.getBoundingClientRect().width, 300);
    } finally {
      pane.remove();
    }
  });

  test('the animal shown follows its study: rebuilt as it changes, cleared for another study (A1, T11)', async () => {
    const n = study();
    let names: StudyInputs | null = n;
    const p = new AnimalPanels({names: () => names, color: () => '#1f77b4', view: treatedView});
    try {
      const shown = p.objectOf('M2');
      grok.shell.o = shown;
      p.refresh();
      expect(grok.shell.o instanceof PmxAnimal && grok.shell.o !== shown, true, 'the same study: built again');
      expect(String(grok.shell.o), 'M2');
      names = study();
      p.refresh();
      expect(grok.shell.o, null, 'another study: cleared, its actions with it');
    } finally {
      grok.shell.o = null;
      p.dispose();
    }
  });

  test('an animal marked from the context panel is shown again there, its actions following the mark', async () => {
    const n = study();
    mark(n.study, [2, 3], 'excluded', 'tumor ulcerated', 'animal M2 (2 rows)');
    const x = new PmxAnimal(n, 'M2');
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = x;
    await delay(100);
    try {
      actionsOf(x).find((a) => a.name === 'Include in the analysis')!.run(x);
      await delay(100);
      const now = grok.shell.o;
      // (the platform may keep the object it shows when the new one is equal: the actions are read from it)
      expect(now instanceof PmxAnimal && (now as PmxAnimal).animal === 'M2', true, String(now));
      expect(names(actionsOf(now as PmxAnimal)).includes('Exclude from the analysis…'), true);
    } finally {
      grok.shell.o = null;
    }
  });

  test('the Actions pane: copying a name in the menu only, a line between marking and the rest', async () => {
    const g = new PmxGroup(study(), 'Vehicle', '#123456');
    const el = withActions(document.createElement('div'), g, GROUP_ACTIONS);
    const links = [...el.querySelectorAll('a, .ui-link')].map((a) => a.textContent?.trim());
    expect(links.includes('Copy the group’s name'), false, links.join('|'));
    expect(links[0], 'Flag every animal of the group…');
    expect(el.querySelectorAll('.pmx-actions-sep').length, 1);
    // the right click still offers it
    expect(names(actionsOf(g)).includes('Copy the group’s name'), true);
  });

  test('a group outside an open study offers no colour; a figure pins and unpins', async () => {
    const g = names(actionsOf(new PmxGroup(study(), 'Vehicle', '#123456')));
    expect(g.includes('Change color…'), false);
    expect(g.includes('Exclude the group from the analysis…'), true);
    let on = false;
    const f = new PmxFigure('The figures', {pin: {on: () => on, set: (x) => on = x}});
    expect(names(actionsOf(f)).join(), 'Pin these figures to compare');
    expect(figurePanel(f).textContent?.includes('Not pinned'), true);
    actionsOf(f)[0].run(f);
    expect(on, true);
    expect(names(actionsOf(f)).join(), 'Unpin: compare with the previous result');
    expect(figurePanel(f).textContent?.includes('New results are drawn beside these'), true);
    // a figure with no chart drawn offers no image
    expect(names(actionsOf(new PmxFigure('growth', {chart: document.createElement('div')}))).length, 0);
  });
});
