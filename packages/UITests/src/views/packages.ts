import {after, before, category} from '@datagrok-libraries/utils/src/test';
// import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


let v: DG.View;

category('Packages', async () => {
  before(async () => {
    v = DG.View.createByType('packages');
  });

  // test('Credentials', async () => {
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
  // });

  after(async () => {
    if (v != null)
      v.close();
  });
}, { owner: 'aparamonov@datagrok.ai' });
