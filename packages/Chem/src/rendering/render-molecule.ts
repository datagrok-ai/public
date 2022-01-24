import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from "openchemlib/full";
import {RDMol} from "../rdkit-api";
import {isMolBlock} from "../utils/chem-utils";
import $ from "cash-dom";
import {convertToRDKit} from "../analysis/r-group-analysis";
import {oclMol} from "../utils/chem-common-ocl";
import {_properties, getRdKitModule, renderer} from "../package";
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import {GridCell, GridCellRenderer} from "datagrok-api/dg";

/** Renders the molecule and returns div with the canvas inside. */
export function renderMolecule(
    molStr: string,
    options?: {renderer?: 'RDKit' | 'OpenChemLib', width?: number, height?: number}): HTMLElement {

  options ??= {};
  options.renderer ??= _properties.Renderer as 'RDKit' | 'OpenChemLib' ?? 'RDKit';
  options.width ??= 200;
  options.height ??= 100;

  let mol: OCL.Molecule | RDMol | null = null;
  let molFile: string;
  let smiles: string;
  isMolBlock(molStr) ? molFile = molStr : smiles = molStr;

  const moleculeHost = ui.canvas(options.width, options.height);
  $(moleculeHost).addClass('chem-canvas');
  // @ts-ignore
  renderer.render(moleculeHost.getContext('2d')!, 0, 0, options.width, options.height, GridCell.fromValue(molStr), null);

  const moreBtn = ui.iconFA(
    'ellipsis-v',
    () => {
      const menu = DG.Menu.popup();
      menu.item('Copy SMILES', () => {
        navigator.clipboard.writeText(smiles);
        grok.shell.info('SMILES copied to clipboard');
      });
      menu.item('Copy Molfile', () => {
        navigator.clipboard.writeText(molFile);
        grok.shell.info('Molfile copied to clipboard');
      });
      menu.item('Sketch', () => {
        const sketcher = new Sketcher();
        isMolBlock(molStr) ? sketcher.setMolFile(molStr) : sketcher.setSmiles(molStr);
        ui.dialog()
          .add(sketcher)
          .show();
      });
      menu.item('Explore', () => {
        grok.shell.o = DG.SemanticValue.fromValueType(molStr, DG.SEMTYPE.MOLECULE)
      });
      menu.show();
    },
    'More',
  );
  $(moreBtn).addClass('chem-mol-view-icon pep-more-icon');

  return ui.divV([moreBtn, moleculeHost], 'chem-mol-box');
}

export function _svgDiv(mol: any) {
  const root = ui.div();
  root.innerHTML = mol.get_svg();
  return root;
}
