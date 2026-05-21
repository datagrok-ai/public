/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ConsensusPharmacophoreApp, loadDemoPdbIds} from './orchestrator';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.app({
    name: 'Consensus Pharmacophore',
    browsePath: 'Bio',
  })
  static async consensusPharmacophoreApp(): Promise<DG.ViewBase> {
    const app = new ConsensusPharmacophoreApp();
    await app.init();
    return app.tableView!;
  }

  @grok.decorators.func({
    name: 'consensusPharmacophoreDemo',
    description: 'Demo: Consensus Pharmacophore on 5 EGFR kinase structures.',
  })
  static async consensusPharmacophoreDemo(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Loading pharmacophore demo...');
    try {
      const app = new ConsensusPharmacophoreApp();
      const pdbIds = await loadDemoPdbIds();
      await app.init(pdbIds);
    } finally { pi.close(); }
  }
}
