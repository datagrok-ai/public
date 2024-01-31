/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {CruddyApp} from "./cruddy";
import {northwindConfig} from "./apps/northwind-app";

export const _package = new DG.Package();

//tags: app
export function northwindDemo() {
  new CruddyApp(northwindConfig).run();
}