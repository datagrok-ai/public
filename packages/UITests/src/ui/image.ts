import {after, before, category, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from './utils';

category('UI: Image', () => {
  let v: DG.View;

  const image = ui.image('https://datagrok.ai/help/visualize/viewers-interaction-main.gif', 400, 200);

  before(async () => {
    v = grok.shell.newView('');
  });

  test('image.root', async () => {
    checkHTMLElement('image', image, v, '.ui-image');
  });

  test('image.click', async () => {
    onClick(image);
  });

  test('image.size', async () => {
    if (image.style.width != '400px' || image.style.height != '200px')
      throw new Error('image size error');
  });

  after(async () => {
    grok.shell.closeAll();
  });

  function onClick(root: HTMLElement): void {
    v.append(root);

    let check = false;
    root.addEventListener('click', function() {
      check = true;
    });

    root.click();

    if (check == false)
      throw new Error(`"${root}": OnClick error`);

    root.remove();
  }
}, {clear: false});
