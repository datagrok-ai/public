import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { caption, enabled, checkHTMLElement } from './utils';

category('UI: Buttons', () => {
  let v: DG.View;
  const buttons = {
    'button': ui.button('', () => { }),
    'bigButton': ui.bigButton('', () => { }),
    'iconButton': ui.button(ui.iconFA(''), () => { }),
  };

  before(async () => {
    v = grok.shell.newView('');
  });

  test('button.root', async () => {
    for (const [key, value] of Object.entries(buttons)) {
      checkHTMLElement(key, value, v, '.ui-btn');
    }
  })

  test('button.click', async () => {
    for (const [key, value] of Object.entries(buttons)) {
      onClick(value);
    }
  })

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });

  function onClick(root: HTMLElement): void {
    v.append(root);

    let check = false;
    root.addEventListener("click", function () {
      check = true;
    });

    root.click();

    if (check == false)
      throw `"${root}": OnClick error`;

    root.remove();
  }

});

