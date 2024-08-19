import {MatchedMolecularPairsViewer} from './mmp-viewer';

declare global {
  interface MatchedMolecularPairsViewer {
    refiltertCliffs(cutoffs: number[], isActiveVar: boolean[], refilter: boolean): void;
  }
}

/**
* Makes the dataset filtering according to user specified cutoffs and 'isActive' checkbox for each activity.
* @param {number[]} cutoffs for each activity the corresponding cuttof is specified by user
* @param {boolean} refilter flag if calculation is required before filtering
*/
MatchedMolecularPairsViewer.prototype.refilterCliffs =
function(cutoffs: number[], isActiveVar: boolean[], refilter: boolean): void {
  if (refilter) {
    for (let i = 0; i < this.cutoffMasks!.length; i ++)
    this.cutoffMasks![i].setAll(false);

    this.totalCutoffMask!.setAll(false);
    this.linesMask!.setAll(false);

    //setting activity associated masks
    for (let i = 0; i < this.lines!.from.length; i++) {
      const activityNumber = this.linesActivityCorrespondance![i];
      const line = this.linesIdxs![i];
      //if 'isActive' checkbox for variable is unset
      //setting line to invisible and continue, otherwise look at cutoff
      if (!isActiveVar[activityNumber]) {
        this.linesMask!.setBit(i, false, false);
        continue;
      }
      if (this.diffs![activityNumber][line] >= cutoffs[activityNumber]) {
        this.cutoffMasks![activityNumber].set(this.lines!.from[i], true, false);
        this.cutoffMasks![activityNumber].set(this.lines!.to[i], true, false);
        this.linesMask!.setBit(i, true, false);
      }
    }

    //setting resulting masks
    //TODO: refine
    for (let j = 0; j < this.totalCutoffMask!.length; j++) {
      let res = false;
      for (let i = 0; i < this.cutoffMasks!.length; i++)
        res = res || this.cutoffMasks![i].get(j);
      this.totalCutoffMask!.set(j, res, false);
    }
  }

  this.parentTable!.filter.copyFrom(this.totalCutoffMask!);
  this.linesRenderer!.linesVisibility = this.linesMask!;
};
