/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {KetcherSketcher} from './ketcher';
export * from './package.g';

export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func({name: 'Ketcher',
    meta: {role: 'moleculeSketcher'},
  })
  static ketcherSketcher() : DG.Widget {
    return new KetcherSketcher();
  }
}
