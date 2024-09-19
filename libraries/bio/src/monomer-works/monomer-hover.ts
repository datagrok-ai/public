import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMonomerHover, ISubstruct, setMonomerHover} from '@datagrok-libraries/chem-meta/src/types';

import {IMonomerLib} from '../types/index';
import {ISeqMonomer} from '../helm/types';
import {HelmTypes} from '../helm/consts';
import {SeqHandler} from '../utils/seq-handler';
import {ALPHABET} from '../utils/macromolecule';
import {helmTypeToPolymerType} from './monomer-works';
import {getMonomersDictFromLib} from './to-atomic-level';
import {monomerSeqToMolfile} from './to-atomic-level-utils';
import {hexToPercentRgb, MonomerHoverLink} from './utils';
import {getMolHighlight} from './seq-to-molfile';

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
        prev.seqColName != seqCol.name || prev.seqPosition != seqMonomer?.position))
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
          getSubstruct: (): ISubstruct | undefined => {
            if (!seqMonomer || seqMonomer.symbol === '*')
              return undefined;

            const seqSH = SeqHandler.forColumn(seqCol);
            const seqSS = seqSH.getSplitted(tableRowIdx);
            const biotype = seqSH.alphabet == ALPHABET.RNA || seqSH.alphabet == ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
            const seqMList: ISeqMonomer[] = wu.count(0).take(seqSS.length)
              .map((posIdx) => { return {position: posIdx, symbol: seqSS.getCanonical(posIdx), biotype: biotype} as ISeqMonomer; })
              .toArray();

            const alphabet = seqSH.alphabet as ALPHABET;
            const polymerType = helmTypeToPolymerType(seqMonomer!.biotype);
            const monomersDict = getMonomersDictFromLib([seqMList], polymerType, alphabet, monomerLib, rdKitModule);
            // Call seq-to-molfile worker core directly
            const molWM = monomerSeqToMolfile(seqMList, monomersDict, alphabet, polymerType);
            const monomerMap = molWM.monomers.get(seqMonomer!.position);
            if (!monomerMap) return {atoms: [], bonds: [], highlightAtomColors: [], highlightBondColors: []};

            const res: ISubstruct = getMolHighlight([monomerMap], monomerLib);
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
  seqCol.temp[MonomerHoverLinksTemp] = mhhList;

  return resLink;
}

export function execMonomerHoverLinks(
  seqGridCell: DG.GridCell, seqMonomer: ISeqMonomer | null
): void {
  const seqCol = seqGridCell.tableColumn!;
  const mhlList = getMonomerHoverLinks(seqCol);
  for (let mhlI = mhlList.length - 1; mhlI >= 0; --mhlI) {
    const mhl = mhlList[mhlI];
    const molGridCol = seqGridCell.grid.col(mhl.targetCol.name);
    if (molGridCol) {
      const handlerRes = mhl.handler(seqGridCell, seqMonomer, molGridCol);
      if (!handlerRes)
        break;
    }
  }
}

export function getMonomerHoverLinks(seqCol: DG.Column<string>): MonomerHoverLink[] {
  return seqCol.temp[MonomerHoverLinksTemp] ?? [];
}
