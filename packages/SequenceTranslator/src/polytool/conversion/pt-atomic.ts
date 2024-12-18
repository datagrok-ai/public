import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export function dealGroups(col: DG.Column<string>): void {
  for (let i = 0; i < col.length; i++) {
    col.set(i, col.get(i)!.replaceAll('undefined', 'H'));
    col.set(i, col.get(i)!.replaceAll('Oh', 'O'));
    col.set(i, col.get(i)!.replaceAll('0.000000 3', '0.000000 0'));
    col.set(i, col.get(i)!.replaceAll('?', 'O'));
    col.set(i, col.get(i)!.replaceAll('0 3\n', '0 0\n'));
    col.set(i, col.get(i)!.replaceAll('RGROUPS=(1 1)', ''));
  }
}

export async function helmToMol(resHelmCol: DG.Column, resList: string[],
  isLinear: boolean[], chiralityEngine: boolean, highlight: boolean, linearize: boolean,
  lib: IMonomerLibBase, rdkit: RDModule, seqHelper: ISeqHelper) {
  const getUnusedName = (df: DG.DataFrame | undefined, colName: string): string => {
    if (!df) return colName;
    return df.columns.getUnusedName(colName);
  };

  const toAtomicLevelRes =
    await seqHelper.helmToAtomicLevel(resHelmCol, chiralityEngine, highlight, lib);

  const resMolCol = toAtomicLevelRes.molCol!;

  const allLinear = isLinear.filter((l) => l).length;
  if (linearize && allLinear > 0) {
    const lin = new Array<string>(allLinear);
    let counter = 0;
    for (let i = 0; i < isLinear.length; i++) {
      if (isLinear[i]) {
        lin[counter] = resList[i];
        counter++;
      }
    }

    const linCol = DG.Column.fromStrings('helm', lin);
    linCol.semType = DG.SEMTYPE.MACROMOLECULE;
    linCol.meta.units = NOTATION.HELM;
    linCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');

    const monomerLibHelper = await getMonomerLibHelper();
    const systemMonomerLib = monomerLibHelper.getMonomerLib();
    try {
      const linear = await _toAtomicLevel(DG.DataFrame.create(0), linCol, systemMonomerLib, seqHelper, rdkit);
      counter = 0;
      for (let i = 0; i < isLinear.length; i++) {
        if (isLinear[i]) {
          resMolCol.set(i, linear!.molCol!.get(counter));
          counter++;
        }
      }
    } catch (e: any) {
      grok.shell.warning('PolyTool was not able to linearize sequences');
    }
  }

  dealGroups(resMolCol);

  return resMolCol;
}
