import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Icons', () => {
  let v: DG.View;

  let icons = ui.icons;
  let iconFA = ui.iconFA('cog');
  let iconSVG = ui.iconSvg('column');
  let iconImage = ui.iconImage('logo', 'http://datagrok.ai/img/logo.svg');

  before(async () => {
    v = grok.shell.newView('');
  });

  test('icons.root', async () => {
    for (const [key, value] of Object.entries(icons)) {
      checkHTMLElement(key, value(() => { }), v, '.grok-icon');
    }
  });

  test('iconFA.root', async () => {
    checkHTMLElement('iconFA', iconFA, v, ['.grok-icon', '.fal']);
  });

  test('iconSVG.root', async () => {
    checkHTMLElement('iconSVG', iconSVG, v, '.svg-icon');
  });

  test('iconImage.root', async () => {
    checkHTMLElement('iconImage', iconImage, v, '.image-icon');
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });

});