import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

import $ from 'cash-dom';

import * as C from '../utils/constants';
import {getAggregatedColumnValues, getStats, Stats} from '../utils/statistics';
import {PeptidesModel} from '../model';
import {getDistributionPanel, getDistributionTable} from '../utils/misc';

export function getDistributionWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  if (!table.selection.anyTrue)
    return new DG.Widget(ui.divText('No distribution'));

  const activityCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY);
  const rowCount = activityCol.length;

  // const setDefaultProperties = (input: DG.InputBase): void => {
  //   input.enabled = !model.isMutationCliffsSelectionEmpty;
  //   $(input.root).find('.ui-input-editor').css('margin', '0px');
  //   $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-left', '5px');
  //   $(input.captionLabel).addClass('ui-label-right');
  // };
  //
  // let defaultValuePos = model.splitByPos;
  // let defaultValueMonomer = model.splitByMonomer;
  // if (!model.isClusterSelectionEmpty && model.isMutationCliffsSelectionEmpty) {
  //   defaultValuePos = false;
  //   defaultValueMonomer = false;
  // }

  const distributionHost = ui.div([], 'd4-flex-wrap');

  const updateDistributionHost = (): void => {
    const res: HTMLDivElement[] = [];
    if (!table.selection.anyTrue)
      res.push(ui.divText('No distribution'));
    else {
      const hist = getActivityDistribution(getDistributionTable(activityCol, model.df.selection, model.getCombinedSelection()));
      const bitArray = BitArray.fromString(table.selection.toBinaryString());
      const mask = DG.BitSet.create(rowCount,
        bitArray.allFalse || bitArray.allTrue ? (_): boolean => true : (i): boolean => bitArray.getBit(i));
      const aggregatedColMap = getAggregatedColumnValues(model.df, model.settings.columns!, {filterDf: true, mask});
      const stats = bitArray.allFalse || bitArray.allTrue ?
        {count: rowCount, pValue: null, meanDifference: 0, ratio: 1, mask: bitArray, mean: activityCol.stats.avg} :
        getStats(activityCol.getRawData(), bitArray);
      const tableMap = getStatsTableMap(stats);
      const resultMap: { [key: string]: any } = {...tableMap, ...aggregatedColMap};
      const distributionRoot = getDistributionPanel(hist, resultMap);
      $(distributionRoot).addClass('d4-flex-col');

      res.push(distributionRoot);
    }

    $(distributionHost).empty().append(res);
  };

  // const splitByPosition = ui.boolInput('Split by position', defaultValuePos, updateDistributionHost);
  // splitByPosition.setTooltip('Constructs distribution for each position separately');
  // setDefaultProperties(splitByPosition);
  // $(splitByPosition.root).css('margin-right', '10px');
  // const splitByMonomer = ui.boolInput('Split by monomer', defaultValueMonomer, updateDistributionHost);
  // splitByMonomer.setTooltip('Constructs distribution for each monomer separately');
  // setDefaultProperties(splitByMonomer);

  // const controlsHost = ui.divH([splitByPosition.root, splitByMonomer.root]);
  // splitByMonomer.fireChanged();
  updateDistributionHost();
  return new DG.Widget(ui.divV([/*controlsHost,*/ distributionHost]));
}

export function getActivityDistribution(table: DG.DataFrame, isTooltip: boolean = false): DG.Viewer<DG.IHistogramLookSettings> {
  const hist = table.plot.histogram({
    filteringEnabled: false,
    valueColumnName: C.COLUMNS_NAMES.ACTIVITY,
    splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
    legendVisibility: 'Never',
    showXAxis: true,
    showColumnSelector: false,
    showRangeSlider: false,
    showBinSelector: false,
    backColor: isTooltip ? '#fdffe5' : '#fffff',
  }) as DG.Viewer<DG.IHistogramLookSettings>;
  hist.root.style.width = 'auto';
  return hist;
}

export function getStatsTableMap(stats: Stats, options: {fractionDigits?: number} = {}): StringDictionary {
  options.fractionDigits ??= 3;
  const tableMap: StringDictionary = {
    'Count': `${stats.count} (${stats.ratio.toFixed(options.fractionDigits)}%)`,
    'Mean difference': stats.meanDifference.toFixed(options.fractionDigits),
    'Mean activity': stats.mean.toFixed(options.fractionDigits),
  };
  if (stats.pValue !== null)
    tableMap['p-value'] = stats.pValue < 0.01 ? '<0.01' : stats.pValue.toFixed(options.fractionDigits);
  return tableMap;
}
