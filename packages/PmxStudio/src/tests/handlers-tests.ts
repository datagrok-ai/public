import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {GROUP_ACTIONS, PmxAnimal, PmxFigure, PmxGroup, actionsOf, figurePanel, registerHandlers, withActions}
  from '../handlers';
import {animalMarks, mark} from '../study-table';
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

  test('an animal: against its twin only where there is one; Exclude, or Include once excluded', async () => {
    const n = study();
    let twin = 0;
    const plain = names(actionsOf(new PmxAnimal(n, 'M1')));
    expect(plain.includes('Show against its twin'), false);
    expect(plain.includes('Exclude from the analysis…') && !plain.includes('Include in the analysis'), true);
    const twinned = new PmxAnimal(n, 'M2', {twin: () => twin++});
    expect(names(actionsOf(twinned))[0], 'Show against its twin');
    actionsOf(twinned)[0].run(twinned);
    expect(twin, 1);
    // excluded: Include is offered, and runs as the old Clear its marks did
    mark(n.study, [2, 3], 'excluded', 'tumor ulcerated', 'animal M2 (2 rows)');
    const after = names(actionsOf(new PmxAnimal(n, 'M2')));
    expect(after.includes('Include in the analysis') && !after.includes('Exclude from the analysis…'), true,
      after.join());
    actionsOf(new PmxAnimal(n, 'M2')).find((a) => a.name === 'Include in the analysis')!.run(new PmxAnimal(n, 'M2'));
    expect(animalMarks(n).get('M2') === undefined, true);
    expect(n.study.col(EXCLUDED)!.get(2), false);
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
