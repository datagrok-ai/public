/* The monomer manager's editor (Bio > Manage > Monomers), a u2 control: its inputs, tabs, R-group rows and
   Save are reached by the library's generic phrases ("R1" input in monomer form, "R-groups" tab); only the
   editor itself and the platform sketcher inside it are named here. */
import {element} from '@datagrok-libraries/bdd';

element('monomer sketcher', {selector: '.monomer-manager-sketcher',
  parts: {'molecule input': 'input[placeholder^="SMILES"]'},
  description: 'the editor\'s sketcher host (DG.chem.Sketcher): its molecule field, its options menu (≡) and the sketcher in it; made with the view\'s editor, with the sketcher of that moment'});
element('monomer form', {selector: '[data-u2-name="monomerEditor"]',
  description: 'the editor of the Manage Monomers view: the library, the sketcher, the monomer\'s inputs on their tabs, an R-groups row per R group and Save, which writes the monomer into the library the view shows'});
