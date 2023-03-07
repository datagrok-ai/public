import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import { structure3d } from './structure3d';
import { identifiers } from './identifiers';

const STRUCT_2D_TAB = '2D';
const STRUCT_3D_TAB = '3D';
const WIDTH = 200;
const HEIGHT = 200;
let lastSelectedTab = STRUCT_2D_TAB;

export function moleculeOverviewWidget(molecule: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  let molblock: string;
  try {
    molblock = _convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }

  let tabControl = ui.tabControl({
    [STRUCT_2D_TAB]: () => {
      return ui.div([
        renderMolecule(molecule, { renderer: 'RDKit', width: WIDTH, height: HEIGHT }), 
        copyAsMolfileIcon(molecule, molblock!)],
        {style: {position: 'relative'}});
    },
    [STRUCT_3D_TAB]: () => {
      return ui.div([ui.wait(async () => await structure3d(molecule, WIDTH, HEIGHT)), copyAsMolfileIcon(molecule, molblock!)], {style: {position: 'relative'}});
    }
  });
  tabControl.currentPane = tabControl.getPane(lastSelectedTab);
  tabControl.root.style.height = `${HEIGHT + 40}px`

  tabControl.onTabChanged.subscribe((_) => {
    lastSelectedTab = tabControl.currentPane.name;
  });

  return new DG.Widget(ui.divV([
    tabControl.root,
    ui.wait(async () => await identifiers(molblock))
  ]));
}

function copyAsMolfileIcon(molecule: string, molblock: string): HTMLButtonElement {
  const copyButton = ui.button(ui.icons.copy(()=>{}, 'Copy as molfile'), ()=>{
    navigator.clipboard.writeText(DG.chem.isMolBlock(molecule) ? molecule : molblock);
  });
  copyButton.classList.add('copy-as-molfile-icon');
  copyButton.style.left = `${WIDTH}px`;
  return copyButton;
}
