/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoView} from './demo-app/demo-app';
// eslint-disable-next-line no-unused-vars
import {expect, expectArray, expectFloat, expectObject} from '@datagrok-libraries/utils/src/test';


export const _package = new DG.Package();


//name: Demo
//tags: app
//description: Interactive demo of major Datagrok capabilities
export function demoApp() {
  grok.shell.addView(new DemoView());
  const pathSegments = window.location.pathname.split('/');
  if (pathSegments.length > 3) {
    const category = pathSegments[3];
    if (category === 'Viewers') {
      const viewerName = pathSegments[4];
      const f = DemoView.findDemoFunc(`${category} | ${viewerName}`);
      f?.apply();
    }
  }
}
