import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import {getHTMLElementbyInnerText} from '../gui/gui-utils';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

let v: DG.View;

category('Packages', async () => {
  before(async () => {
    v = DG.View.createByType('packages');
  });

  test('Changelog', async () => {
    grok.shell.addView(v);
    await delay(1000);
    const packages = v.root.getElementsByClassName('grok-package-name');
    const uiTestsPackage = Array.from(packages).find((p) => (p as HTMLElement).textContent!.includes('UI Tests'));
    if (uiTestsPackage == null)
      throw new Error('UI Tests package not found in the view');
    const version = uiTestsPackage.getElementsByClassName('grok-package-version')[0].textContent;
    const event = new KeyboardEvent('keydown', {
        key: 'F4',
        code: 'F4',
        keyCode: 115,
        which: 115,
        altKey: false,
        ctrlKey: false,
        shiftKey: false,
        metaKey: false,
        bubbles: true,
        cancelable: true,
    });
    document.dispatchEvent(event);
    (uiTestsPackage as HTMLElement).click();
    await delay(1000);
    const changelogDiv = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Release Notes');
    if (changelogDiv == null)
      throw new Error('Release Notes section not found');
    const changelogText = changelogDiv.textContent;
    if (changelogText == null)
      throw new Error('Release Notes section contains no text');
    expect(changelogText.includes(version!), true);
  }, {timeout: 5000});

  //test('Credentials', async () => {
  //   grok.shell.addView(v);
  //   const divs = v.root.querySelectorAll('div.grok-package-name');
  //   divs.forEach((d) => {
  //     if (d.firstChild!.textContent == 'API Tests') {
  //       console.log(d);
  //       const e = d.ownerDocument.createEvent('MouseEvents');
  //       e.initMouseEvent('contextmenu', true, true,
  //         element.ownerDocument.defaultView, 1, 0, 0, 0, 0, false,
  //         false, false, false,2, null);


  //       return !element.dispatchEvent(e);
  //     }
  //   });
  //});

  after(async () => {
    if (v != null)
      v.close();
  });
});
