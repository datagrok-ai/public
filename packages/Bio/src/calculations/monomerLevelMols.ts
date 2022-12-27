import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from '../utils/constants';
import {getHelmMonomers} from '../package';
import {getSplitter, getStats} from '@datagrok-libraries/bio';

const V2000_ATOM_NAME_POS = 31;

export async function getMonomericMols(mcol: DG.Column,
  pattern: boolean = false, monomersDict?: Map<string, string>): Promise<DG.Column> {
  const separator: string = mcol.tags[C.TAGS.SEPARATOR];
  const units: string = mcol.tags[DG.TAGS.UNITS];
  const splitter = getSplitter(units, separator);
  let molV3000Array;
  monomersDict ??= new Map();
  const monomers = units === 'helm' ?
    getHelmMonomers(mcol) : Object.keys(getStats(mcol, 0, splitter).freq).filter((it) => it !== '');

  for (let i = 0; i < monomers.length; i++) {
    if (!monomersDict.has(monomers[i]))
      monomersDict.set(monomers[i], `${monomersDict.size + 1}`);
  }

  if (units === 'helm') {
    molV3000Array = await grok.functions.call('HELM:getMolFiles', {col: mcol});
    molV3000Array = changeV2000ToV3000(molV3000Array, monomersDict, pattern);
  } else {
    molV3000Array = new Array<string>(mcol.length);
    for (let i = 0; i < mcol.length; i++) {
      const sequenceMonomers = splitter(mcol.get(i)).filter((it) => it !== '');
      const molV3000 = molV3000FromNonHelmSequence(sequenceMonomers, monomersDict, pattern);
      molV3000Array[i] = molV3000;
    }
  }
  return DG.Column.fromStrings('monomericMols', molV3000Array);
}

function molV3000FromNonHelmSequence(
  monomers: Array<string>, monomersDict: Map<string, string>, pattern: boolean = false) {
  let molV3000 = `
  Datagrok macromolecule handler

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
`;

  molV3000 += `M  V30 COUNTS ${monomers.length} ${monomers.length ? monomers.length - 1 : 0} 0 0 0\n`;
  molV3000 += 'M  V30 BEGIN ATOM\n';

  for (let atomRowI = 0; atomRowI < monomers.length; atomRowI++) {
    molV3000 += pattern ?
      `M  V30 ${atomRowI + 1} R${monomersDict.get(monomers[atomRowI])} 0.000 0.000 0 0\n` :
      `M  V30 ${atomRowI + 1} At 0.000 0.000 0 0 MASS=${monomersDict.get(monomers[atomRowI])}\n`;
  }

  molV3000 += 'M  V30 END ATOM\n';
  molV3000 += 'M  V30 BEGIN BOND\n';

  for (let bondRowI = 0; bondRowI < monomers.length - 1; bondRowI++)
    molV3000 += `M  V30 ${bondRowI + 1} 1 ${bondRowI + 1} ${bondRowI + 2}\n`;

  molV3000 += 'M  V30 END BOND\n';
  molV3000 += 'M  V30 END CTAB\n';
  molV3000 += 'M  END';
  return molV3000;
}

function changeV2000ToV3000(mols: DG.Column, dict: Map<string, string>, pattern: boolean = false): Array<string> {
  const molsArray = new Array<string>(mols.length);
  for (let i = 0; i < mols.length; i++) {
    let curPos = 0;
    let endPos = 0;
    let molV3000 = `
    Datagrok macromolecule handler

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
`;

    const mol = mols.get(i);
    curPos = mol.indexOf('\n', curPos) + 1;
    curPos = mol.indexOf('\n', curPos) + 1;
    curPos = mol.indexOf('\n', curPos) + 1;

    const atomMonomerCounts = parseInt(mol.substring(curPos, curPos + 3));
    const bondMonomerCounts = parseInt(mol.substring(curPos + 3, curPos + 6));

    molV3000 += `M  V30 COUNTS ${atomMonomerCounts} ${bondMonomerCounts} 0 0 0\n`;
    molV3000 += 'M  V30 BEGIN ATOM\n';

    for (let atomRowI = 0; atomRowI < atomMonomerCounts; atomRowI++) {
      curPos = mol.indexOf('\n', curPos) + 1 + V2000_ATOM_NAME_POS;
      endPos = mol.indexOf(' ', curPos);
      const monomerName: string = mol.substring(curPos, endPos);
      molV3000 += pattern ?
        `M  V30 ${atomRowI + 1} R${dict.get(monomerName)} 0.000 0.000 0 0\n` :
        `M  V30 ${atomRowI + 1} At 0.000 0.000 0 0 MASS=${dict.get(monomerName)}\n`;
    }

    molV3000 += 'M  V30 END ATOM\n';
    molV3000 += 'M  V30 BEGIN BOND\n';

    for (let bondRowI = 0; bondRowI < bondMonomerCounts; bondRowI++) {
      curPos = mol.indexOf('\n', curPos) + 1;
      const firstMonomer = parseInt(mol.substring(curPos, curPos + 3).trim());
      const secondMonomer = parseInt(mol.substring(curPos + 3, curPos + 6).trim());
      const order = parseInt(mol.substring(curPos + 6, curPos + 9).trim());

      molV3000 += `M  V30 ${bondRowI + 1} ${order} ${firstMonomer} ${secondMonomer}\n`;
    }

    molV3000 += 'M  V30 END BOND\n';
    molV3000 += 'M  V30 END CTAB\n';
    molV3000 += 'M  END';
    molsArray[i] = molV3000;
  }

  return molsArray;
}
