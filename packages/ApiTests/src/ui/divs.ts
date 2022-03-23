import {after, before, category, delay, expect, test} from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from "./utils";

category('UI - div', () => {
  const htmlElement = ui.divText('test');
  let v: DG.View;
  const divs: {[key: string]: {'element': HTMLDivElement, 'selectors': string[]}} = {
    'block': {'element': ui.block([htmlElement]), 'selectors': ['.ui-div', '.ui-block']},
    'block25': {'element': ui.block25([htmlElement]), 'selectors': ['.ui-div', '.ui-block', '.ui-block-25']},
    'block50': {'element': ui.block50([htmlElement]), 'selectors': ['.ui-div', '.ui-block', '.ui-block-50']},
    'block75': {'element': ui.block75([htmlElement]), 'selectors': ['.ui-div', '.ui-block', '.ui-block-75']},
    'box': {'element': ui.box(htmlElement), 'selectors': ['.ui-box']},
    'boxFixed': {'element': ui.boxFixed(htmlElement), 'selectors': ['.ui-box', '.ui-box-fixed']},
    'card': {'element': ui.card(htmlElement), 'selectors': ['.d4-item-card', '.ui-div']},
    'div': {'element': ui.div([htmlElement]), 'selectors': ['.ui-div']},
    'divH': {'element': ui.divH([htmlElement, htmlElement]), 'selectors': ['.d4-flex-row', '.ui-div']}, 
    'divText': {'element': ui.divText('test'), 'selectors': ['div']},
    'divV': {'element': ui.divV([htmlElement, htmlElement]), 'selectors': ['.d4-flex-col', '.ui-div']},
    'info': {'element': ui.info(htmlElement), 'selectors': ['.grok-info-bar-container', '.ui-div']},
    'panel': {'element': ui.panel([htmlElement]), 'selectors': ['.ui-div', '.ui-panel']},
    'splitH': {'element': ui.splitH([htmlElement, htmlElement]), 'selectors': ['.ui-box', '.ui-split-h']},
    'splitV': {'element': ui.splitV([htmlElement, htmlElement]), 'selectors': ['.ui-box', '.ui-split-v']},
    'buttonsInput': {'element': ui.buttonsInput(), 'selectors': ['.ui-input-root', '.ui-input-buttons']},
    'form': {'element': ui.form(), 'selectors': ['.ui-form']},
    'image': {'element': ui.image('https://dev.datagrok.ai/favicon.png', 42, 42), 'selectors': ['.ui-image']},
    'inputs': {'element': ui.inputs([ui.floatInput('test', 42)]), 'selectors': ['.ui-form']},
    'narrowForm': {'element': ui.narrowForm(), 'selectors': ['.ui-form', '.ui-form-condensed']},
    'loader': {'element': ui.loader(), 'selectors': ['.grok-loader']},
    'wait': {'element': ui.wait(async ()=> ui.div()), 'selectors': ['.grok-wait']},
    'waitBox': {'element': ui.waitBox(async ()=> ui.div()), 'selectors': ['.grok-wait', '.ui-box']},
  };

  before(async () => {
    v = grok.shell.newView();
  });

  test('div.root', async () => {
    for (const [key, value] of Object.entries(divs)) {
      checkHTMLElement(key, value['element'], v, value['selectors']);
    }
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
