import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import $ from 'cash-dom';

import * as C from '../utils/constants';
import {getStats, Stats} from '../utils/statistics';
import {PeptidesModel} from '../model';
import {getStatsSummary, prepareTableForHistogram} from '../utils/misc';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

const allConst = 'All';
const otherConst = 'Other';

export function getDistributionWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  const activityCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const activityColData = activityCol.getRawData();
  const rowCount = activityCol.length;
  const selectionObject = model.invariantMapSelection;
  const clustersColName = model.settings.clustersColumnName;
  let clustersProcessedObject: string[] = [];
  if (clustersColName)
    clustersProcessedObject = Object.values(model.clusterSelection).flat();

  const positions = Object.keys(selectionObject);
  let monomerStr = allConst;
  let otherStr = '';

  const updateDistributionHost = (): void => {
    model.splitByPos = splitByPosition.value!;
    model.splitByMonomer = splitByMonomer.value!;
    const res: HTMLDivElement[] = [];
    if (splitByPosition.value && splitByMonomer.value) {
      otherStr = otherConst;
      for (const position of positions) {
        const monomerList = selectionObject[position];
        if (monomerList.length === 0)
          continue;

        const posCol = table.getCol(position);
        const posColCategories = posCol.categories;
        const posColData = posCol.getRawData();

        for (const monomer of monomerList) {
          const labels = getDistributionLegend(`${position} : ${monomer}`, otherStr);

          const monomerCategoryIndex = posColCategories.indexOf(monomer);
          const mask = DG.BitSet.create(rowCount, (i) => posColData[i] === monomerCategoryIndex);
          const distributionTable = DG.DataFrame.fromColumns(
            [activityCol, DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask)]);
          const hist = getActivityDistribution(prepareTableForHistogram(distributionTable));

          const stats = model.monomerPositionStats[position]![monomer]!;
          const tableMap = getStatsTableMap(stats);

          const aggregatedColMap = model.getAggregatedColumnValues({filterDf: true, mask});

          const resultMap = {...tableMap, ...aggregatedColMap};
          const distributionRoot = getStatsSummary(labels, hist, resultMap);
          $(distributionRoot).addClass('d4-flex-col');

          res.push(distributionRoot);
        }
      }
    } else if (splitByPosition.value) {
      otherStr = otherConst;
      for (const position of positions) {
        const monomerList = selectionObject[position];
        if (monomerList.length === 0)
          continue;

        monomerStr = `${position}: {${monomerList.join(', ')}}`;
        const labels = getDistributionLegend(monomerStr, otherStr);

        const posCol = table.getCol(position);
        const posColCategories = posCol.categories;
        const posColData = posCol.getRawData();
        const monomerIndexesList = monomerList.map((monomer) => posColCategories.indexOf(monomer));
        const mask = DG.BitSet.create(rowCount, (i) => monomerIndexesList.includes(posColData[i]));
        const splitCol = DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask);

        const aggregatedColMap = model.getAggregatedColumnValues({filterDf: true, mask});

        const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
        const hist = getActivityDistribution(prepareTableForHistogram(distributionTable));

        const bitArray = BitArray.fromUint32Array(rowCount, splitCol.getRawData() as Uint32Array);
        const stats = getStats(activityColData, bitArray);
        const tableMap = getStatsTableMap(stats);

        const resultMap = {...tableMap, ...aggregatedColMap};
        const distributionRoot = getStatsSummary(labels, hist, resultMap);
        $(distributionRoot).addClass('d4-flex-col');

        res.push(distributionRoot);
      }
    } else if (splitByMonomer.value) {
      const reversedSelectionObject: {[monomer: string]: string[]} = {};
      const monomers = [];
      for (const position of positions) {
        for (const monomer of selectionObject[position]) {
          if (!reversedSelectionObject.hasOwnProperty(monomer)) {
            reversedSelectionObject[monomer] = [position];
            monomers.push(monomer);
            continue;
          }
          if (!reversedSelectionObject[monomer].includes(position))
            reversedSelectionObject[monomer].push(position);
        }
      }

      otherStr = otherConst;
      for (const monomer of monomers) {
        const posList = reversedSelectionObject[monomer];
        const posColList = posList.map((pos) => table.getCol(pos));
        const posColCategoriesList = posColList.map((posCol) => posCol.categories);
        const posColDataList = posColList.map((posCol) => posCol.getRawData());
        const monomerCategoryIndexList = posColCategoriesList.map((posColCategories) => posColCategories.indexOf(monomer));

        monomerStr = `${monomer}: {${posList.join(', ')}}`;
        const labels = getDistributionLegend(monomerStr, otherStr);

        const mask = DG.BitSet.create(rowCount,
          (i) => posColDataList.some((posColData, j) => posColData[i] === monomerCategoryIndexList[j]));
        const aggregatedColMap = model.getAggregatedColumnValues({filterDf: true, mask});

        const splitCol = DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask);
        const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
        const hist = getActivityDistribution(prepareTableForHistogram(distributionTable));

        const bitArray = BitArray.fromUint32Array(rowCount, splitCol.getRawData() as Uint32Array);
        const stats = getStats(activityColData, bitArray);
        const tableMap = getStatsTableMap(stats);

        const resultMap: {[key: string]: any} = {...tableMap, ...aggregatedColMap};
        const distributionRoot = getStatsSummary(labels, hist, resultMap);
        $(distributionRoot).addClass('d4-flex-col');

        res.push(distributionRoot);
      }
    } else {
      const splitCol = table.col(C.COLUMNS_NAMES.SPLIT_COL);
      if (!splitCol)
        res.push(ui.divText('No distribution'));
      else {
        otherStr = '';
        if (Object.values(selectionObject).some((selectedAar) => selectedAar.length !== 0) ||
          clustersProcessedObject.length !== 0) {
          monomerStr = '';
          for (const position of positions) {
            const monomerList = selectionObject[position];
            if (monomerList.length !== 0)
              monomerStr += `${position}: {${monomerList.join(', ')}}; `;
          }
          if (clustersProcessedObject.length !== 0)
            monomerStr += `Clusters: ${clustersProcessedObject.join(', ')}`;
          otherStr = otherConst;
        }
        const labels = getDistributionLegend(monomerStr, otherStr);
        const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
        const hist = getActivityDistribution(prepareTableForHistogram(distributionTable));
        const bitArray = BitArray.fromUint32Array(rowCount, splitCol.getRawData() as Uint32Array);
        const mask = DG.BitSet.create(rowCount,
          bitArray.allFalse ? (_): boolean => true : (i): boolean => bitArray.getBit(i));
        const aggregatedColMap = model.getAggregatedColumnValues({filterDf: true, mask});
        const stats = bitArray.allFalse ? {count: rowCount, pValue: null, meanDifference: 0, ratio: 1, mask: bitArray} :
          getStats(activityColData, bitArray);
        const tableMap = getStatsTableMap(stats);
        const resultMap: {[key: string]: any} = {...tableMap, ...aggregatedColMap};
        const distributionRoot = getStatsSummary(labels, hist, resultMap);
        $(distributionRoot).addClass('d4-flex-col');

        res.push(distributionRoot);
      }
    }
    $(distributionHost).empty().append(res);
  };

  const setDefaultProperties = (input: DG.InputBase): void => {
    input.enabled = !model.isMonomerPositionSelectionEmpty;
    $(input.root).find('.ui-input-editor').css('margin', '0px');
    $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-left', '5px');
  };

  let defaultValuePos = model.splitByPos;
  let defaultValueMonomer = model.splitByMonomer;
  if (!model.isClusterSelectionEmpty && model.isMonomerPositionSelectionEmpty) {
    defaultValuePos = false;
    defaultValueMonomer = false;
  }

  const splitByPosition = ui.boolInput('', defaultValuePos, updateDistributionHost);
  splitByPosition.addPostfix('Split by position');
  splitByPosition.setTooltip('Constructs distribution for each position separately');
  setDefaultProperties(splitByPosition);
  $(splitByPosition.root).css('margin-right', '10px');
  const splitByMonomer = ui.boolInput('', defaultValueMonomer, updateDistributionHost);
  splitByMonomer.addPostfix('Split by monomer');
  splitByMonomer.setTooltip('Constructs distribution for each monomer separately');
  setDefaultProperties(splitByMonomer);

  const controlsHost = ui.divH([splitByPosition.root, splitByMonomer.root]);
  const distributionHost = ui.div([], 'd4-flex-wrap');
  splitByMonomer.fireChanged();

  return new DG.Widget(ui.divV([controlsHost, distributionHost]));
}

export function getActivityDistribution(table: DG.DataFrame, isTooltip: boolean = false,
): DG.Viewer<DG.IHistogramLookSettings> {
  const hist = table.plot.histogram({
    filteringEnabled: false,
    valueColumnName: C.COLUMNS_NAMES.ACTIVITY_SCALED,
    splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
    legendVisibility: 'Never',
    showXAxis: true,
    showColumnSelector: false,
    showRangeSlider: false,
    showBinSelector: !isTooltip,
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
  };
  if (stats.pValue !== null)
    tableMap['p-value'] = stats.pValue < 0.01 ? '<0.01' : stats.pValue.toFixed(options.fractionDigits);
  return tableMap;
}

export function getDistributionLegend(thisLabel: string, otherLabel: string = ''): HTMLDivElement {
  return ui.divV([
    ui.label(thisLabel, {style: {color: DG.Color.toHtml(otherLabel.length === 0 ? DG.Color.blue : DG.Color.orange)}}),
    ui.label(otherLabel, {style: {color: DG.Color.toHtml(DG.Color.blue)}})]);
}
