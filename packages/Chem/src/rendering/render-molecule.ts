import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from "openchemlib/full";
import {RDMol} from "../rdkit-api";
import {isMolBlock} from "../utils/chem-utils";
import $ from "cash-dom";
import {convertToRDKit} from "../analysis/r-group-analysis";
import {oclMol} from "../utils/chem-common-ocl";
import {_properties, getRdKitModule} from "../package";
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;


/** Renders the molecule and returns div with the canvas inside. */
export function renderMolecule(
    molStr: string,
    options: {renderer?: 'RDKit' | 'OpenChemLib', width?: number, height?: number}): HTMLElement {

  options.renderer ??= _properties.Renderer as 'RDKit' | 'OpenChemLib';
  options.width ??= 200;
  options.height ??= 150;

  let mol: OCL.Molecule | RDMol | null = null;
  let molFile: string;
  let smiles: string;
  isMolBlock(molStr) ? molFile = molStr : smiles = molStr;

  const moleculeHost = ui.canvas(options.width, options.height);
  $(moleculeHost).addClass('chem-canvas');

  switch (options.renderer) {
    case 'RDKit':
      try {
        mol = getRdKitModule().get_mol(convertToRDKit(molStr)!);
        (mol as RDMol).draw_to_canvas(moleculeHost, options.width, options.height);
        molFile ??= (mol as RDMol).get_molblock();
        smiles ??= (mol as RDMol).get_smiles();
      } finally {
        (mol as RDMol)?.delete();
      }
      break;
    case 'OpenChemLib':
      mol = oclMol(molStr);
      OCL.StructureView.drawMolecule(moleculeHost, mol as OCL.Molecule);
      molFile ??= mol.toMolfile();
      smiles ??= mol.toSmiles();
      break;
    default:
      throw new Error(`Renderer '${options.renderer}' is not supported.`);
  }

  const moreBtn = ui.iconFA(
    'ellipsis-v',
    () => {
      const menu = DG.Menu.popup();
      menu.item('Copy SMILES', () => {
        navigator.clipboard.writeText(smiles);
        grok.shell.info('SMILES copied!');
      });
      menu.item('Copy Molfile', () => {
        navigator.clipboard.writeText(molFile);
        grok.shell.info('Molfile copied!');
      });
      menu.item('Sketch', () => {
        const sketcher = new Sketcher();
        isMolBlock(molStr) ? sketcher.setMolFile(molStr) : sketcher.setSmiles(molStr);
        ui.dialog()
          .add(sketcher)
          .show();
      });
      menu.item('Explore', () => {
        grok.shell.o = DG.SemanticValue.fromValueType(molStr, 'Molecule')
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
