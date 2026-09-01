/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IMonomerLibProvider} from '@datagrok-libraries/bio/src/types/monomer-library';
import {DomainDbLibraryProvider} from './domain-library-provider';

export * from './package.g';

export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func({
    name: 'Monomers Domain DB Provider',
    outputs: [{type: 'object', name: 'result'}],
    meta: {role: 'monomer-lib-provider'},
  })
  static async getMonomerDomainDbProvider(): Promise<IMonomerLibProvider> {
    return DomainDbLibraryProvider.getInstance();
  }
}
