/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {addSubstructProvider, getMonomerHover, ISubstruct, setMonomerHover} from '@datagrok-libraries/chem-meta/src/types';

import {IMonomerLibBase} from '../types/monomer-library';
import {ISeqMonomer} from '../helm/types';
import {MonomerHoverLink} from './utils';
import {getMolHighlight} from './seq-to-molfile';
import {ISeqHelper} from '../utils/seq-helper';
import {AsmLayoutName, getSeqSrcLayout, MonomerMapProvider} from './assembler';

export const MonomerHoverLinksTemp = 'MonomerHoverLinks';

export function addMonomerHoverLink(seqColTemp: any, resLink: MonomerHoverLink) {
  let mhhList = seqColTemp[MonomerHoverLinksTemp];
  if (!mhhList)
    mhhList = seqColTemp[MonomerHoverLinksTemp] = [];
  mhhList.push(resLink);
  seqColTemp[MonomerHoverLinksTemp] = mhhList;
}

/** Invalidates the grids that may show {@link molCol}: the table views of its data frame and the last grid a
 * hover came from. */
function invalidateMolGrids(molCol: DG.Column<string>, lastGrid: DG.Grid | null): void {
  lastGrid?.invalidate();
  const df = molCol.dataFrame;
  if (!df) return;
  for (const tv of grok.shell.tableViews)
    if (tv.dataFrame?.id === df.id) tv.grid?.invalidate();
}

/**
 * Builds a monomer hover link between sequence and molecule columns. The monomer maps come from the macromolecule
 * assembler per row, on demand (D6): the first hover of a row requests its map and redraws the molecule cell when it
 * arrives (no highlight until then); with the `.sequence-src-highlight-monomers` tag on, the rows the renderer draws
 * get their maps requested and their grids redrawn when they arrive. The maps are built with the layout in the
 * molecule column's `.sequence-src-layout` tag (any value but 'linear' is non-linear). Molfiles not made by the
 * assembler (stored by the old JavaScript conversion) do not match its maps (D7).
 * @param {DG.Column<string>} seqCol the sequence column the molecules were made from
 * @param {DG.Column<string>} molCol the molecule column
 * @param {IMonomerLibBase} monomerLib the library the molecules were made with
 * @param {ISeqHelper} seqHelper
 * @param {RDModule} rdKitModule deprecated and unused: the maps come from the assembler
 * @param {boolean} throughPOM the layout for a molecule column without the layout tag: non-linear when true
 * @return {Promise<MonomerHoverLink>} the registered link
 */
export async function buildMonomerHoverLink(
  seqCol: DG.Column<string>, molCol: DG.Column<string>,
  monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, rdKitModule?: RDModule | null, throughPOM: boolean = false,
): Promise<MonomerHoverLink> {
  const defaultLayout: AsmLayoutName = throughPOM ? 'nonlinear' : 'linear';
  const maps = new MonomerMapProvider(seqCol, monomerLib, seqHelper, () => getSeqSrcLayout(molCol, defaultLayout));

  let lastGrid: DG.Grid | null = null;
  let invalidateTimer: number | null = null;
  /** Maps of the whole-molecule highlight arrive one by one; one redraw for those arriving together */
  const scheduleInvalidate = (): void => {
    if (invalidateTimer !== null) return;
    invalidateTimer = window.setTimeout(() => {
      invalidateTimer = null;
      invalidateMolGrids(molCol, lastGrid);
    }, 50);
  };

  const resLink: MonomerHoverLink = {
    targetCol: molCol,
    handler: (seqGridCell: DG.GridCell, seqMonomer: ISeqMonomer | null, targetGridCol: DG.GridColumn): boolean => {
      const grid = targetGridCol.grid;
      lastGrid = grid;
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
        const seqPosition = seqMonomer.position;
        setMonomerHover({
          gridCell: targetGridCell,
          dataFrameId: seqCol.dataFrame.id,
          gridRowIdx: gridRowIdx,
          seqColName: seqCol.name,
          // The assembler keys maps by the sequence position the renderers report (a FASTA/separator nucleotide's
          // sugar, base and phosphate share the base's position), so no remapping is needed
          seqPosition: seqPosition ?? -1,
          getSubstruct: (): ISubstruct | undefined => { // Gets monomer highlight
            if (!seqMonomer || seqMonomer.symbol === '*')
              return undefined;

            const molMonomerMap = maps.get(tableRowIdx); // undefined until the map arrives
            if (!molMonomerMap)
              return undefined;

            const monomerMap = molMonomerMap.get(seqPosition); // single monomer
            if (!monomerMap) return {atoms: [], bonds: [], highlightAtomColors: [], highlightBondColors: []};

            const res: ISubstruct = getMolHighlight([monomerMap], monomerLib);
            return res;
          }
        });

        // Request the hovered row's map first; redraw the molecule cell when it arrives, if still hovered
        maps.get(tableRowIdx, () => {
          const cur = getMonomerHover();
          if (cur && cur.dataFrameId == seqCol.dataFrame?.id && cur.gridRowIdx == gridRowIdx &&
            cur.seqColName == seqCol.name && cur.seqPosition == seqPosition)
            targetGridCell.render();
        }, true);

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

      // Prefetch: the rows the renderer draws get their maps requested; their grids redraw when the maps arrive
      const molMonomerMap = maps.get(tableRowIdx, scheduleInvalidate);
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
