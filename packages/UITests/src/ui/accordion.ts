import {after, before, category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from './utils';

category('UI: Accordion', () => {
  let v: DG.View;
  let acc: DG.Accordion ;

  before(async () => {
    acc = ui.accordion()
    v = grok.shell.newView('');
    v.append(acc.root);
    acc.addPane('New pane',()=> ui.div(''))
    acc.addPane('New pane3',()=> ui.div(''))
  });

  test('accordion.root', async () => {
    checkHTMLElement('accordion', acc.root, v, '.d4-accordion');
  });

  test('accordion.addPane', async () => {
    acc.addPane('New pane1', () => ui.div(), false);
    checkHTMLElement('accordion pane', acc.root, v, '.d4-accordion-pane');
  });

  test('accordion.getPane', async () => {
    if (acc.getPane('New pane') == undefined)
      throw new Error('getPane error');
  });

  test('pane.name', async () => {
    if (acc.panes[0].name != 'New pane')
      throw new Error('pane.name error');
  });

  test('pane.root', async () => {
    checkHTMLElement('Accordion Pane root', acc.panes[0].root, v, '.d4-accordion-pane');
  });

  test('pane.expanded', async () => {
    if (acc.panes[0].expanded != false)
      throw new Error('pane.name error');
  });

  test('accordion.removePane', async () => {
    debugger
    acc.removePane(acc.getPane('New pane3'));
    if (acc.getPane('New pane3') != undefined)
      throw new Error('getPane error');
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {clear: false});
