/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {northwindConfig} from "./apps/northwind-app";

import {CruddyApp} from "./cruddy/app";
import {chemblConfig} from "./apps/chembl-app";

export const _package = new DG.Package();

//tags: app
//meta.browsePath: Dev
export function northwindDemo() {
  new CruddyApp(northwindConfig).run();
}

//tags: app
export function chemblDemo() {
  new CruddyApp(chemblConfig).run();
}