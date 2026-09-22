/* The datasets Helm's features open (the package's own samples on the stand, published with it)
   and the HELM editor: a full-screen dialog with no title (`name="dialog-"`), whose own parts come
   from the hwe library and carry `data-testid`, not the u2 contract — so they are named here, each
   scoped to the editor. Platform names (grid, context panel, ...) are reserved. */
import {dataset, element} from '@datagrok-libraries/bdd';

dataset('helm-showcase', {path: 'System:AppData/Helm/samples/helm-showcase.csv', aliases: ['helm showcase'],
  description: '53 curated HELM cases: Name, Category, HELM — row 1 PEPTIDE1{A.C}$$$$, row 2 the 10-mer A..L'});
dataset('HELM', {path: 'System:AppData/Helm/samples/HELM.csv', aliases: ['helm sample'],
  description: 'the package\'s large HELM peptide sample: HELM and Activity'});

const EDITOR = 'HELM editor';
const testid = (id: string) => `[data-testid="${id}"]`;
element(EDITOR, {selector: `.d4-dialog-full-screen:has(${testid('app-root')})`, aliases: ['HELM web editor'],
  description: 'the full-screen HELM editor dialog a HELM cell opens (double-click, or Edit Helm... on its context menu)'});
const inEditor: Record<string, [string, string[]?]> = {
  'editor canvas': [testid('editor-svg'), ['drawing area']],
  'drawn monomers': ['[data-testid^="canvas-atom-"]', ['drawn monomer']],
  'Sequence tab': [testid('tab-sequence')],
  'HELM tab': [testid('tab-helm')],
  'Properties tab': [testid('tab-properties')],
  'notation pane': [testid('notation-pane-content'), ['HELM notation']],
  'notation error': [testid('notation-pane-error')],
  'formula field': [testid('properties-formula')],
  'molecular weight field': [testid('properties-mw')],
  'extinction coefficient field': [testid('properties-extinction')],
  'palette search': [testid('palette-search')],
  'Peptides palette tab': [testid('palette-tab-PEPTIDE')],
  'RNA palette tab': [testid('palette-tab-RNA')],
  'Favorites palette tab': [testid('palette-tab-Favorites')],
  'favorites empty note': [testid('palette-favorites-empty')],
  'monomer tiles': ['[data-testid^="palette-tile-"]', ['monomer tile']],
  'G monomer tile': [testid('palette-tile-G')],
  'Aca monomer tile': [testid('palette-tile-Aca')],
  'RNA builder': [testid('rna-builder')],
  'triplets': ['[data-testid^="palette-triplet-"]', ['triplet']],
  'editor status': [testid('status-mode'), ['editor status bar']],
  'undo button': [testid('toolbar-undo')],
  'redo button': [testid('toolbar-redo')],
  'clean layout button': [testid('toolbar-clean')],
};
for (const [name, [selector, aliases]] of Object.entries(inEditor))
  element(name, {selector, aliases, in: EDITOR});
