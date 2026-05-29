/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {loadDemoPdbIds} from './orchestrator';
import {WizardShell} from './wizard-shell';
import {StepId} from './wizard-state';

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
    // Wizard-based entry: builds a custom DG.ViewBase with a top step rail +
    // persistent Mol* viewer. Replaces the old TableView-with-toolbox-buttons UI.
    const shell = new WizardShell();
    return shell.getView();
  }

  @grok.decorators.func({
    name: 'consensusPharmacophoreDemo',
    description: 'Demo: Consensus Pharmacophore on 5 EGFR kinase structures.',
  })
  static async consensusPharmacophoreDemo(): Promise<void> {
    // Demo entry — load the EGFR PDB list into the wizard and auto-run the
    // full pipeline so the user sees the final consensus pharmacophore
    // without clicking through each step.
    const pi = DG.TaskBarProgressIndicator.create('Loading pharmacophore demo...');
    try {
      const pdbIds = await loadDemoPdbIds();
      const shell = new WizardShell(pdbIds);
      grok.shell.addView(shell.getView());
      // Run on a microtask so the view is fully attached to the DOM before
      // the preview begins manipulating Mol*. Without this, the Mol* viewer
      // mounts into an off-screen container on some Chrome / DG combinations.
      setTimeout(() => {
        shell.runFullPipeline().catch((e: any) => {
          grok.shell.warning(`Demo pipeline failed: ${e?.message ?? e}. ` +
            'You can still drive the app manually.');
        });
      }, 0);
    } finally { pi.close(); }
  }
}

// Re-export StepId so external scripts can reference step IDs symbolically.
export {StepId};
