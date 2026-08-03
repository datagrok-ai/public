/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';
import {LRUCache} from 'lru-cache';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {addSubstructProvider, getMonomerHover, ISubstruct, setMonomerHover} from '@datagrok-libraries/chem-meta/src/types';

import {IMonomerLibBase} from '../types/monomer-library';
import {HelmType, ISeqMonomer, PolymerType} from '../helm/types';
import {HelmTypes, PolymerTypes} from '../helm/consts';
import {ALPHABET} from '../utils/macromolecule';
import {buildRolesForHelmRna, getMonomersDictFromLib} from './to-atomic-level';
import {monomerSeqToMolfile} from './to-atomic-level-utils';
import {MonomerHoverLink} from './utils';
import {getMolHighlight} from './seq-to-molfile';
import {MonomerMap} from './types';
import {IHelmToMolfileConverter, ISeqHelper} from '../utils/seq-helper';

export const MonomerHoverLinksTemp = 'MonomerHoverLinks';

export function addMonomerHoverLink(seqColTemp: any, resLink: MonomerHoverLink) {
  let mhhList = seqColTemp[MonomerHoverLinksTemp];
  if (!mhhList)
    mhhList = seqColTemp[MonomerHoverLinksTemp] = [];
  mhhList.push(resLink);
  seqColTemp[MonomerHoverLinksTemp] = mhhList;
}

/**
 * Builds a monomer hover link between sequence and molecule columns
 * throughPOM specifies if the sequence is translated through seqhelper using POM (default is false)
 */
export async function buildMonomerHoverLink(
  seqCol: DG.Column<string>, molCol: DG.Column<string>,
  monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, rdKitModule: RDModule, throughPOM: boolean = false
): Promise<MonomerHoverLink> {
  const seqSH = seqHelper.getSeqHandler(seqCol);
  const isNucleotide = (seqSH.alphabet == ALPHABET.RNA || seqSH.alphabet == ALPHABET.DNA);
  const isNucleotideHelmSequence = seqSH.isHelm() && isNucleotide;
  // if the conversion goes through polymer object model, it is aware that DNA/RNA need sugars and phosphates, therefore, for fasta
  // of length 15, it will return 45 monomers. so we need to handle this case separately
  const isNucleotideNonHelmSequence = !seqSH.isHelm() && isNucleotide;

  const getSeqMonomerCorrectedPosition = (pos?: number) => {
    if (pos == undefined) return null;
    // HELM RNA is split into per-position sugar/base/phosphate monomers, so the
    // renderer reports a flat monomer index and the monomer map is keyed by that
    // same flat index — both for the POM converter and (since the map is now
    // built with roles) for the linear path. No division needed in either case.
    if (isNucleotideHelmSequence)
      return pos;
    // Non-HELM DNA/RNA: the linear map keys by base position; the POM converter
    // expands each base into a sugar/base/phosphate triple, so map the base to
    // its base entry (the middle of its triple).
    if (isNucleotideNonHelmSequence)
      return throughPOM ? pos * 3 + 1 : pos;
    return pos;
  };

  const helmToMolfileConverter: IHelmToMolfileConverter | null = throughPOM ? await seqHelper.getHelmToMolfileConverter(monomerLib) : null;
  function buildMonomerMap(seqCol: DG.Column<string>, tableRowIdx: number): MonomerMap {
    const seqSH = seqHelper.getSeqHandler(seqCol);
    if (!throughPOM) {
      const seqSS = seqSH.getSplitted(tableRowIdx);
      const gi = seqSS.graphInfo;

      // Polymer type: for HELM read it from the parsed graph rather than the
      // alphabet — modified RNA can mis-detect as ALPHABET.UN. The whole row is
      // assumed to be the first chain's type. Non-HELM falls back to the alphabet.
      let polymerType: PolymerType;
      if (seqSH.isHelm() && gi?.polymerTypes && gi.polymerTypes.length > 0)
        polymerType = gi.polymerTypes[0];
      else {
        const a = seqSH.alphabet as ALPHABET;
        polymerType = (a === ALPHABET.RNA || a === ALPHABET.DNA) ? PolymerTypes.RNA : PolymerTypes.PEPTIDE;
      }

      // Coerce a UN alphabet to RNA for the RNA polymer path (only used for the
      // bases-only default sugar/phosphate selection; irrelevant in triples mode).
      let alphabet = seqSH.alphabet as ALPHABET;
      if (polymerType === PolymerTypes.RNA && alphabet !== ALPHABET.RNA && alphabet !== ALPHABET.DNA)
        alphabet = ALPHABET.RNA;

      // HELM RNA is assembled from per-position sugar/base/phosphate triples;
      // the monomer map must be built with the same roles (and multi-chain
      // stacking) the real conversion uses, otherwise every monomer is treated
      // as a base and the atom map does not match the rendered molfile.
      const isHelmRna = seqSH.isHelm() && polymerType === PolymerTypes.RNA;
      const biotype: HelmType = isHelmRna ? (HelmTypes.NUCLEOTIDE as HelmType) : seqSH.defaultBiotype;
      const seqMList: ISeqMonomer[] = wu.count(0).take(seqSS.length)
        .map((posIdx) => { return {position: posIdx, symbol: seqSS.getCanonical(posIdx), biotype: biotype} as ISeqMonomer; })
        .toArray();

      const chainStarts = (gi?.disjointSeqStarts && gi.disjointSeqStarts.length > 0) ?
        gi.disjointSeqStarts.slice() : [0];
      const roles = isHelmRna ?
        buildRolesForHelmRna([seqMList], monomerLib, polymerType, [chainStarts])[0] : undefined;

      const monomersDict = getMonomersDictFromLib([seqMList], [roles], polymerType, alphabet, monomerLib, rdKitModule);
      // Call seq-to-molfile worker core directly
      const molWM = monomerSeqToMolfile(seqMList, monomersDict, alphabet, polymerType, roles, chainStarts);
      return molWM.monomers;
    } else {
      const helm = seqSH.getHelm(tableRowIdx);
      const molWM = seqHelper.helmToAtomicLevelSingle(helm, helmToMolfileConverter!, false, false); // skip mol beautification step. not needed.
      return molWM.monomers;
    }
  }

  const monomerMapLruCache = new LRUCache<string, MonomerMap>({max: 100});

  function getMonomerMap(seqCol: DG.Column<string>, tableRowIdx: number): MonomerMap | null {
    const seq = seqCol.get(tableRowIdx);
    if (seq == null) return null;

    let resMonomerMap = monomerMapLruCache.get(seq);
    if (!resMonomerMap)
      monomerMapLruCache.set(seq, resMonomerMap = buildMonomerMap(seqCol, tableRowIdx));

    return resMonomerMap;
  }

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
          //prev.gridCell.grid?.invalidate();
          prev.gridCell.render();
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
          seqPosition: getSeqMonomerCorrectedPosition(seqMonomer?.position) ?? -1,
          getSubstruct: (): ISubstruct | undefined => { // Gets monomer highlight
            if (!seqMonomer || seqMonomer.symbol === '*')
              return undefined;

            const molMonomerMap = getMonomerMap(seqCol, tableRowIdx);
            if (!molMonomerMap)
              return undefined;

            const monomerMap = molMonomerMap.get(getSeqMonomerCorrectedPosition(seqMonomer?.position)!); // single monomer
            if (!monomerMap) return {atoms: [], bonds: [], highlightAtomColors: [], highlightBondColors: []};

            const res: ISubstruct = getMolHighlight([monomerMap], monomerLib);
            return res;
          }
        });

        // TODO: Invalidate targetGridCell
        //grid.invalidate();
        targetGridCell.render();
      }

      return true;
    },
    /* ISubstructProvider.*/getSubstruct: (tableRowIdx: number | null,): ISubstruct | undefined => { // Gets whole molecule highlight
      if (molCol.getTag(ChemTags.SEQUENCE_SRC_HL_MONOMERS) != 'true') return undefined;
      if (tableRowIdx == null) return undefined;
      const seq = seqCol.get(tableRowIdx);
      if (!seq) return undefined;

      const molMonomerMap = getMonomerMap(seqCol, tableRowIdx);
      if (!molMonomerMap) return undefined;
      const res: ISubstruct = getMolHighlight(molMonomerMap.values(), monomerLib);
      return res;
    }
  };

  addMonomerHoverLink(seqCol.temp, resLink);
  addSubstructProvider(molCol.temp, resLink);

  return resLink;
}

export function execMonomerHoverLinks(
  seqGridCell: DG.GridCell, seqMonomer: ISeqMonomer | null
): void {
  try {
    const seqCol = seqGridCell.tableColumn;
    if (!seqCol) return;

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
  } catch (e) {
    console.error(e);
  }
}

export function getMonomerHoverLinks(seqCol: DG.Column<string>): MonomerHoverLink[] {
  return seqCol.temp[MonomerHoverLinksTemp] ?? [];
}
