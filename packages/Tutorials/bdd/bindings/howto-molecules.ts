/* Demo tables the Molecules how-to clips need besides the library's (spgi, smiles, mol1K): Chem's own demo files
   for the analyses that want an activity column. Opened off camera, like the library's datasets. */
import {dataset, kind} from '@datagrok-libraries/bdd';

// as Chem's bdd project names it (Chem/bdd/bindings/elements.ts), which the Tutorials project does not load
kind('calculator', {
  selector: '.biochem-calc-nav-item',
  match: ['text'],
  editorSelector: 'input[type="checkbox"]',
  description: 'a calculator row of the Chemical Properties dialog, by its name ("Chemical Properties", "logP")',
});

dataset('sar_small', {path: 'System:DemoFiles/chem/sar_small.csv', aliases: ['sar small'],
  description: '200 molecules of one series (smiles) with assay columns: In-vivo Activity, Assay #1(IC50), cLogP...'});
dataset('activity_cliffs', {path: 'System:DemoFiles/chem/activity_cliffs.csv', aliases: ['activity cliffs demo'],
  description: '269 kinase inhibitors (smiles) with an Activity column, plus a precomputed embedding and cliff flags'});
