import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMonomerHover, ISubstruct, setMonomerHover} from '@datagrok-libraries/chem-meta/src/types';

import {IMonomerLib} from '../types/index';
import {ISeqMonomer} from '../helm/types';
import {SeqHandler} from '../utils/seq-handler';
import {ALPHABET} from '../utils/macromolecule';
import {helmTypeToPolymerType} from './monomer-works';
import {getMonomersDictFromLib} from './to-atomic-level';
import {monomerSeqToMolfile} from './to-atomic-level-utils';
import {getGridColByTableCol} from '../utils/grid';
import {hexToPercentRgb, MonomerHoverLink} from './utils';

export const MonomerHoverLinksTemp = 'MonomerHoverLinks';

export function buildMonomerHoverLink(
  seqCol: DG.Column<string>, molCol: DG.Column<string>, monomerLib: IMonomerLib, rdKitModule: RDModule
): MonomerHoverLink {
  const resLink: MonomerHoverLink = {
    targetCol: molCol,
    handler: (seqGridCell: DG.GridCell, seqMonomer: ISeqMonomer | null, targetGridCol: DG.GridColumn): boolean => {
      const grid = targetGridCol.grid;
      const tableRowIdx = seqGridCell.tableRowIndex!;
      const gridRowIdx = seqGridCell.gridRow;
      const targetGridCell = grid.cell(targetGridCol.name, gridRowIdx);

      const prev = getMonomerHover();
      if (!prev || (prev && (prev.dataFrameId != seqCol.dataFrame.id || prev.gridRowIdx != gridRowIdx ||
        prev.seqColName != seqCol.name || prev.seqPosition != seqMonomer?.position ||
        prev.molColName != molCol.name))
      ) {
        if (prev) {
          setMonomerHover(null);
          prev.gridCell.grid?.invalidate();
        }
        if (!seqMonomer) {
          setMonomerHover(null);
          return true;
        }

        setMonomerHover({
          gridCell: targetGridCell,
          dataFrameId: seqCol.dataFrame.id,
          gridRowIdx: gridRowIdx,
          seqColName: seqCol.name,
          seqPosition: seqMonomer ? seqMonomer.position : -1,
          molColName: molCol.name,
          getSubstruct: (): ISubstruct | undefined => {
            if (!seqMonomer || seqMonomer.symbol === '*')
              return undefined;

            const seqSH = SeqHandler.forColumn(seqCol);
            const seqSS = seqSH.getSplitted(tableRowIdx);
            const seqCMList = wu.count(0).take(seqSS.length).map((posIdx) => seqSS.getCanonical(posIdx)).toArray();

            const alphabet = seqSH.alphabet as ALPHABET;
            const polymerType = helmTypeToPolymerType(seqMonomer!.biotype);
            const monomersDict = getMonomersDictFromLib([seqCMList], polymerType, alphabet, monomerLib, rdKitModule);
            // Call seq-to-molfile worker core directly
            const molWM = monomerSeqToMolfile(seqCMList, monomersDict, alphabet, polymerType);
            const monomerMap = molWM.monomers.get(seqMonomer!.position);
            if (!monomerMap) return {atoms: [], bonds: [], highlightAtomColors: [], highlightBondColors: []};

            const hlAtoms: { [key: number]: number[] } = {};
            const hlBonds: { [key: number]: number[] } = {};
            const wem = monomerLib.getWebEditorMonomer(seqMonomer.biotype, seqMonomer.symbol)!;
            const mColorStr = wem.backgroundcolor!;
            const wemColorA = hexToPercentRgb(mColorStr ?? DG.Color.mouseOverRows) ?? [1.0, 0.0, 0.0, 0.7];
            for (const mAtom of monomerMap.atoms) hlAtoms[mAtom] = wemColorA;
            for (const mBond of monomerMap.bonds) hlBonds[mBond] = wemColorA;

            const res: ISubstruct = {
              atoms: monomerMap.atoms,
              bonds: monomerMap.bonds,
              highlightAtomColors: hlAtoms,
              highlightBondColors: hlBonds,
            };
            return res;
          }
        });

        // TODO: Invalidate targetGridCell
        grid.invalidate();
      }

      return true;
    },
  };

  let mhhList = seqCol.temp[MonomerHoverLinksTemp];
  if (!mhhList)
    mhhList = seqCol.temp[MonomerHoverLinksTemp] = [];
  mhhList.push(resLink);

  return resLink;
}

export function execMonomerHoverLinks(
  seqGridCell: DG.GridCell, seqMonomer: ISeqMonomer | null
): void {
  const seqCol = seqGridCell.tableColumn!;
  const mhlList = getMonomerHoverLinks(seqCol);
  for (let mhlI = mhlList.length - 1; mhlI >= 0; --mhlI) {
    const mhl = mhlList[mhlI];
    const molGridCol = getGridColByTableCol(seqGridCell.grid, mhl.targetCol);
    if (molGridCol)
      if (!mhl.handler(seqGridCell, seqMonomer, molGridCol)) break;
  }
}

export function getMonomerHoverLinks(seqCol: DG.Column<string>): MonomerHoverLink[] {
  return seqCol.temp[MonomerHoverLinksTemp] ?? [];
}
