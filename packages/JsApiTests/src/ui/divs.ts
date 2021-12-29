import {after, before, category, delay, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from "./utils";

category('UI - div', () => {
  const htmlElement = ui.divText('test');
  let v: DG.View;
  const divs = {
    'block': ui.block([htmlElement]),
    'block25': ui.block25([htmlElement]),
    'block50': ui.block50([htmlElement]),
    'block75': ui.block75([htmlElement]),
    'box': ui.box(htmlElement),
    'boxFixed': ui.boxFixed(htmlElement),
    'card': ui.card(htmlElement),
    'div': ui.div([htmlElement]),
    'divH': ui.divH([htmlElement, htmlElement]),
    'divText': ui.divText('test'),
    'divV': ui.divV([htmlElement, htmlElement]),
    'info': ui.info(htmlElement),
    'panel': ui.panel([htmlElement]),
    'splitH': ui.splitH([htmlElement, htmlElement]),
    'splitV': ui.splitV([htmlElement, htmlElement]),
  };

  before(async () => {
    v = grok.shell.newView();
  });

  test('div.root', async () => {
    for (const [key, value] of Object.entries(divs)) {
      checkHTMLElement(key, value, v, '.ui-div');
    }
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
