/* The monomer manager's editor (Bio > Manage > Monomers), a u2 control: its inputs, the R-group menus,
   Details toggle and Save are reached by the library's generic phrases ("Symbol" input in monomer form, "Br"
   menu item); the editor, the platform sketcher inside it and the R-group chips are named here. */
import {element, kind} from '@datagrok-libraries/bdd';

element('monomer sketcher', {selector: '.monomer-manager-sketcher',
  parts: {'molecule input': 'input[placeholder^="SMILES"]'},
  description: 'the editor\'s sketcher host (DG.chem.Sketcher): its molecule field, its options menu (≡) and the sketcher in it; made with the view\'s editor, with the sketcher of that moment'});
element('monomer form', {selector: '[data-u2-name="monomerEditor"]',
  description: 'the editor of the Manage Monomers view: the library, the sketcher, the monomer\'s symbol, name and types two to a row, a chip per R group, and, pinned at the bottom, Save, the line saying what it does or why it cannot yet, and the Details toggle, whose natural analog, ID, colors and properties open above it'});
kind('R-group chip', {selector: '.monomer-rgroup-chip', match: ['item'],
  description: 'an R group of the monomer in the editor, by its label ("R2"), reading "R2 OH": invalid (red) while it cannot be saved, the reason in its tooltip; a click opens its menu of the common caps ("Br" menu item, the current one checked) and "Edit definition…", which opens the fields a library stores ("Cap SMILES" input) under the chips; a second click closes the menu'});
