import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectBoolGetSet, expectNoThrow} from '../helpers';

category('AI: Widgets: Accordion JS API', () => {
  test('create + addPane + panes + getPane + expanded round-trip', async () => {
    const acc = DG.Accordion.create();
    expect(acc instanceof DG.Accordion, true);
    const pane = acc.addPane('alpha', () => ui.divText('alpha-content'));
    expect(pane instanceof DG.AccordionPane, true);
    expect(pane.name, 'alpha');
    expect(acc.panes.length, 1);
    const fetched = acc.getPane('alpha');
    expect(fetched instanceof DG.AccordionPane, true);
    expect(fetched.name, 'alpha');
    pane.expanded = true;
    expect(pane.expanded, true);
    pane.expanded = false;
    expect(pane.expanded, false);
  });

  test('addTitle + addCountPane build structure', async () => {
    const acc = ui.accordion();
    expect(acc instanceof DG.Accordion, true);
    expectNoThrow(() => acc.addTitle(ui.divText('Section title')));
    const counted = acc.addCountPane('beta', () => ui.divText('beta-content'), () => 7);
    expect(counted instanceof DG.AccordionPane, true);
    expect(counted.name, 'beta');
    expect(acc.panes.length, 1);
  });

  test('autoHideTabHeader get/set round-trip', async () => {
    expectBoolGetSet(DG.Accordion.create(), 'autoHideTabHeader');
  });

  test('context get/set round-trip', async () => {
    const acc = DG.Accordion.create();
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', [1, 2, 3])]);
    acc.context = df;
    expect(acc.context instanceof DG.DataFrame, true);
    expect((acc.context as DG.DataFrame).rowCount, 3);
  });

  test('removePane shrinks panes; end finalizes', async () => {
    const acc = DG.Accordion.create();
    const p1 = acc.addPane('one', () => ui.divText('1'));
    const p2 = acc.addPane('two', () => ui.divText('2'));
    expect(acc.panes.length, 2);
    acc.removePane(p1);
    expect(acc.panes.length, 1);
    expect(acc.getPane('two') instanceof DG.AccordionPane, true);
    expectNoThrow(() => acc.end());
    expect(p2.name, 'two');
  });

  test('boundary: empty accordion + double-end', async () => {
    const acc = DG.Accordion.create();
    expect(acc.panes.length, 0);
    expectNoThrow(() => acc.end());
    expectNoThrow(() => acc.end());
    expect(acc.panes.length, 0);
  });
}, {owner: 'agolovko@datagrok.ai'});
