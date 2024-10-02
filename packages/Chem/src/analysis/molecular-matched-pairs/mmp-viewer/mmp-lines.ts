import {ILineSeries} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {PaletteCodes} from './palette';
import {MMPA} from '../mmp-analysis/mmpa';


export function createLines(mmpa: MMPA, palette: PaletteCodes) : {
    linesIdxs: Uint32Array,
    lines: ILineSeries,
    linesActivityCorrespondance: Uint32Array
  } {
  let allCount = 0;
  for (let i = 0; i < mmpa.initData.activitiesCount; i++)
    allCount += mmpa.allCasesBased.activityPairsIdxs[i].trueCount();

  const pointsFrom = new Uint32Array(allCount);
  const pointsTo = new Uint32Array(allCount);
  const linesIdxs = new Uint32Array(allCount);
  const colors = new Array<string>(allCount);
  const linesActivityCorrespondance = new Uint32Array(allCount);

  let pairsCounter = 0;
  for (let j = 0; j < mmpa.initData.activitiesCount; j++) {
    for (let i = -1; (i = mmpa.allCasesBased.activityPairsIdxs[j].findNext(i)) !== -1;) {
      pointsFrom[pairsCounter] = mmpa.allCasesBased.molNumFrom[i];
      pointsTo[pairsCounter] = mmpa.allCasesBased.molNumTo[i];
      linesIdxs[pairsCounter] = i;
      colors[pairsCounter] = palette.rgbCut[j];
      linesActivityCorrespondance[pairsCounter] = j;
      pairsCounter++;
    }
  }

  const lines: ILineSeries = {
    from: pointsFrom,
    to: pointsTo,
    drawArrows: true,
    colors: colors,
    arrowSize: 10,
    width: 0.5,
  };

  return {linesIdxs, lines, linesActivityCorrespondance};
}
