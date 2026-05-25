/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ConsensusPharmacophoreApp, loadDemoPdbIds} from './orchestrator';
import {DEFAULT_OPTIONS} from './orchestrator-types';

// Bundle the package stylesheet — webpack's style-loader injects it into
// <head> at runtime so all `cp-*` selectors used by orchestrator.ts and
// cluster-picker-dialog.ts actually render. See css/consensus-pharmacophore.css
// for the full rule set.
import '../css/consensus-pharmacophore.css';

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
    // Demo entry — the user clicks once and ends up at the consensus
    // pharmacophore rendering, with no further interaction needed. Previously
    // this only preloaded the textarea and dumped the user on a half-loaded
    // screen; the new flow runs the full pipeline so the demo is what it
    // claims to be. We use the Consensus *preview* path (not runPipeline) on
    // purpose: the demo PDBs are all active-state EGFR kinases, so the
    // cluster picker would auto-skip anyway, and the preview path doesn't
    // open the picker dialog at all — keeps the demo zero-click.
    const pi = DG.TaskBarProgressIndicator.create('Loading pharmacophore demo...');
    try {
      const app = new ConsensusPharmacophoreApp();
      const pdbIds = await loadDemoPdbIds();
      await app.init(pdbIds);
      // Run on a microtask so the TableView is fully attached to the DOM
      // before the preview begins manipulating Mol*. Without this, Mol*
      // mounts into an off-screen container on some Chrome / DG combinations.
      setTimeout(() => {
        app.previewConsensus(pdbIds, DEFAULT_OPTIONS).catch((e) => {
          grok.shell.warning(`Demo pipeline failed: ${e?.message ?? e}. ` +
            'You can still drive the app manually.');
        });
      }, 0);
    } finally { pi.close(); }
  }
}
