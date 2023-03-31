import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {_properties, renderer, convertMolNotation} from '../package';
import {RDKitCellRenderer} from './rdkit-cell-renderer';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import { _convertMolNotation } from '../utils/convert-notation-utils';
import { isSmarts } from '../utils/mol-creation_rdkit';

/** Renders the molecule and returns div with the canvas inside. */
export function renderMolecule(
  molStr: string,
  options?: {renderer?: 'RDKit' | 'OpenChemLib',
  width?: number, height?: number, popupMenu?: boolean}): HTMLElement {
  options ??= {};
  options.renderer ??= _properties.Renderer as 'RDKit' | 'OpenChemLib' ?? 'RDKit';
  options.width ??= 200;
  options.height ??= 100;
  options.popupMenu ??= true;

  const moleculeHost = ui.canvas(options.width, options.height);

  $(moleculeHost).addClass('chem-canvas');
  const r = window.devicePixelRatio;
  moleculeHost.width = options.width * r;
  moleculeHost.height = options.height * r;
  moleculeHost.style.width = (options.width).toString() + 'px';
  moleculeHost.style.height = (options.height).toString() + 'px';

  const renderFunctions = DG.Func.find({meta: {chemRendererName: options.renderer}});
  if (renderFunctions.length > 0) {
    renderFunctions[0].apply().then((rendndererObj) => {
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, options!.width, options!.height,
        DG.GridCell.fromValue(molStr));
    });
  }

  const div = ui.divV([], 'chem-mol-box');

  if (options.popupMenu) {
    const moreBtn = ui.iconFA(
      'ellipsis-v',
      () => {
        const menu = DG.Menu.popup();
        menu.item('Copy SMILES', () => {
          const smiles = !DG.chem.isMolBlock(molStr) && !isSmarts(molStr) ? molStr :
            _convertMolNotation(molStr, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, getRdKitModule());
          navigator.clipboard.writeText(smiles);
          grok.shell.info('SMILES copied to clipboard');
        });
        menu.item('Copy Molfile', () => {
          const molFile = DG.chem.isMolBlock(molStr) ? molStr : _convertMolNotation(molStr, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, getRdKitModule());
          navigator.clipboard.writeText(molFile);
          grok.shell.info('Molfile copied to clipboard');
        });
        menu.item('Sketch', () => {
          const sketcher = new DG.chem.Sketcher();
          DG.chem.isMolBlock(molStr) ? sketcher.setMolFile(molStr) : sketcher.setSmiles(molStr);
          ui.dialog()
            .add(sketcher)
            .show();
        });
        menu.item('Explore', () => {
          grok.shell.o = DG.SemanticValue.fromValueType(molStr, DG.SEMTYPE.MOLECULE);
        });
        menu.show();
      },
      'More',
    );
    $(moreBtn).addClass('chem-mol-view-icon pep-more-icon');
    div.append(moreBtn);
  }
  div.append(moleculeHost);

  return div;
}

export function _svgDiv(mol: any): HTMLDivElement {
  const root = ui.div();
  root.innerHTML = mol.get_svg();
  return root;
}
