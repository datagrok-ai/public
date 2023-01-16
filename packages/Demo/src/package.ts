/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DemoView} from "./demo-app/demo-app";


export let _package = new DG.Package();


//name: Demo
//tags: app
//description: Interactive demo of major Datagrok capabilities
export function demoApp() {
  grok.shell.addView(new DemoView());
}