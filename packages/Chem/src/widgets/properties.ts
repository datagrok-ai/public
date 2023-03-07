import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol} from '../utils/chem-common-ocl';
import {div} from "datagrok-api/ui";
import $ from 'cash-dom';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import { MOL_FORMAT } from '../constants';

export function getMoleculeCharge(mol: OCL.Molecule) {
  const atomsNumber = mol.getAllAtoms();
  let moleculeCharge = 0;
  for (let atomIndx = 0; atomIndx <= atomsNumber; ++atomIndx) {
    moleculeCharge += mol.getAtomCharge(atomIndx);
  }
  return moleculeCharge;
}

export function propertiesWidget(semValue: DG.SemanticValue<string>): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    semValue.value = _convertMolNotation(semValue.value, 'unknown', MOL_FORMAT.SMILES, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  let host = div();
  try {
    var mol = oclMol(semValue.value);
  } catch {
    return new DG.Widget(ui.divText('Could not analyze properties'));
  }

  function prop(name: string, type: DG.ColumnType, extract: (mol: OCL.Molecule) => any) {
    var addColumnIcon = ui.iconFA('plus', () => {
      let molCol: DG.Column<string> = semValue.cell.column;
      semValue.cell.dataFrame.columns
        .addNew(semValue.cell.dataFrame.columns.getUnusedName(name), type)
        .init((i) => {
          try {
            if (molCol.isNone(i)) return null;
            const mol = oclMol(molCol.get(i)!);
            return extract(mol);
          }
          catch (_) {
            return null;
          }
        });
    }, `Calculate ${name} for the whole table`);

    ui.tools.setHoverVisibility(host, [addColumnIcon]);
    $(addColumnIcon)
      .css('color', '#2083d5')
      .css('position', 'absolute')
      .css('top', '2px')
      .css('left', '-12px')
      .css('margin-right', '5px');

    return ui.divH([addColumnIcon, extract(mol)], { style: {'position': 'relative'}});
  }

  let map = {
    'Formula': prop('Formula', DG.TYPE.STRING, (m) => m.getMolecularFormula().formula),
    'MW': prop('MW', DG.TYPE.FLOAT, (m) => m.getMolecularFormula().absoluteWeight),
    'HBA': prop('HBA', DG.TYPE.INT, (m) => new OCL.MoleculeProperties(m).acceptorCount),
    'HBD': prop('HBD', DG.TYPE.INT, (m) => new OCL.MoleculeProperties(m).donorCount),
    'LogP': prop('LogP', DG.TYPE.FLOAT, (m) => new OCL.MoleculeProperties(m).logP),
    'LogS': prop('LogS', DG.TYPE.FLOAT, (m) => new OCL.MoleculeProperties(m).logS),
    'PSA': prop('PSA', DG.TYPE.FLOAT, (m) => new OCL.MoleculeProperties(m).polarSurfaceArea),
    'Rotatable bonds': prop('Rotatable bonds', DG.TYPE.INT, (m) => new OCL.MoleculeProperties(m).rotatableBondCount),
    'Stereo centers': prop('Stereo centers', DG.TYPE.INT, (m) => new OCL.MoleculeProperties(m).stereoCenterCount),
    'Molecule charge': prop('Molecule charge', DG.TYPE.INT, (m) => getMoleculeCharge(m)),
  };

  host.appendChild(ui.tableFromMap(map));
  return new DG.Widget(host);
}