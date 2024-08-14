import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {MMP_NAMES} from '../mmp-constants';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MatchedMolecularPairsViewer} from './mmp-viewer';

declare global {
  interface MatchedMolecularPairsViewer {
    refreshPair(dkitModule: RDModule): Promise<void>;
    refreshFilterAllFragments(): void;
    refilterAllFragments(rowChanged: boolean): void;
    pinPair(idx: number): void;
  }
}

/**
 * Prepares all the entities to show for selected pair.
 * @param {RDModule} rdkitModule RDkit module instance
 */
MatchedMolecularPairsViewer.prototype.refreshPair =
async function(rdkitModule: RDModule) : Promise<void> {
  const progressBarPairs = DG.TaskBarProgressIndicator.create(`Refreshing pairs...`);

  this.transPairsMask!.setAll(false);
  const diffFromSubstrCol = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.STRUCT_DIFF_FROM_NAME);
  const diffToSubstrCol = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.STRUCT_DIFF_TO_NAME);
  const diffFrom = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.FROM);
  const diffTo = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.TO);

  const [idxPairs, cases] = this.findSpecificRule(diffFromSubstrCol);
  await this.recoverHighlights(cases, diffFrom, diffTo, diffFromSubstrCol, diffToSubstrCol, rdkitModule);

  this.transPairsGrid!.dataFrame.filter.copyFrom(this.transPairsMask!);

  this.transPairsGrid!.setOptions({
    pinnedRowColumnNames: [],
    pinnedRowValues: [],
  });
  if (idxPairs >= 0) {
    this.transPairsGrid!.setOptions({
      pinnedRowValues: [idxPairs.toString()],
      pinnedRowColumnNames: [MMP_NAMES.PAIRNUM],
    });
  }

  progressBarPairs.close();
};


MatchedMolecularPairsViewer.prototype.refreshFilterAllFragments = function(): void {
  // const consistsBitSet: DG.BitSet = DG.BitSet.create(this.transFragmentsGrid!.dataFrame.rowCount);
  // consistsBitSet.setAll(true);
  this.transFragmentsGrid!.dataFrame.filter.copyFrom(this.fragmentsMask!);
};

/**
* Gets visible all fragment pairs according to selected molecule in the table.
* @param {boolean} rowChanged flag if row was changed
*/
MatchedMolecularPairsViewer.prototype.refilterAllFragments =
function(rowChanged: boolean): void {
  let idxTrue = -1;
  if (rowChanged) {
    const idx = this.parentTable!.currentRowIdx;
    this.transFragmentsMask!.setAll(false);

    for (let i = 0; i < this.mmpRules!.rules.length; i++) {
      for (let j = 0; j < this.mmpRules!.rules[i].pairs.length; j++) {
        const fs = this.mmpRules!.rules[i].pairs[j].firstStructure;
        if (idx == fs) {
          if (idxTrue == -1)
            idxTrue = i;
          this.transFragmentsMask!.set(i, true, false);
          break;
        }
      }
    }
  }

  if (this.enableFilters) {
    this.transFragmentsGrid!.dataFrame.filter.copyFrom(this.transFragmentsMask!);
    this.transFragmentsGrid!.invalidate();
    if (rowChanged)
      this.transFragmentsGrid!.table.currentRowIdx = idxTrue;
  }
};

/**
* pin the pair from "Pairs" dataframe to parent table
* @param {number} idx index of pair
*/
MatchedMolecularPairsViewer.prototype.pinPair = function(idx: number): void {
  const columns = this.transPairsGrid?.dataFrame.columns;
  const idxFrom: number = columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(idx);
  const idxToTo: number = columns!.byName(MMP_NAMES.PAIRNUM_TO).get(idx);
  const molFrom = this.parentTable!.columns.byName(this.molecules!).get(idxFrom);
  const molTo = this.parentTable!.columns.byName(this.molecules!).get(idxToTo);
  const grid = grok.shell.tv.grid;
  const indexesPairs = this.transPairsMask?.getSelectedIndexes();
  const indexesAllFrom = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(ip));
  const indexesAllTo = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_TO).get(ip));
  indexesAllFrom?.forEach((i) => grid.dataFrame.selection.set(i, true));
  indexesAllTo?.forEach((i) => grid.dataFrame.selection.set(i, true));

  grid.setOptions({
    pinnedRowValues: [molFrom, molTo],
    pinnedRowColumnNames: [this.molecules!, this.molecules!],
  });
};
