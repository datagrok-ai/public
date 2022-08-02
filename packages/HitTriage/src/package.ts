/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp, HitTriageTemplate} from "./hit-triage-app";

export let _package = new DG.Package();


//tags: app
//name: Hit Triage
export async function hitTriageApp() {
  new HitTriageApp(HitTriageTemplate.demo());
}