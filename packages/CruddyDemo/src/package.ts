/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {northwindConfig} from "./apps/northwind-app";

import {CruddyApp} from "./cruddy/app";
import {chemblConfig} from "./apps/chembl-app";
export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions{
  @grok.decorators.app({
      'browsePath': 'Dev' 
  })
  static northwindDemo() {
  
    new CruddyApp(northwindConfig).run();
  }


  @grok.decorators.app({})
  static chemblDemo() {
  
    new CruddyApp(chemblConfig).run();
  }
}