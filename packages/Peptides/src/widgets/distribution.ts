import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import $ from 'cash-dom';

import * as C from '../utils/constants';
import {getStats, Stats} from '../utils/statistics';
import {PeptidesModel} from '../model';

const allConst = 'All';
const otherConst = 'Other';

export function getDistributionWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  const activityScaledCol = table.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
  const rowCount = activityScaledCol.length;
  const selectionObject = model.mutationCliffsSelection;
  const clustersObject = model.logoSummarySelection;
  const positions = Object.keys(selectionObject);
  const positionsLen = positions.length;
  let aarStr = allConst;
  let otherStr = '';
  const useSelectedStr = model.isPeptideSpaceChangingBitset;

  const updateDistributionHost = (): void => {
    model.splitByPos = splitByPosition.value!;
    model.splitByAAR = splitByAAR.value!;
    const res: HTMLDivElement[] = [];
    if (splitByPosition.value && splitByAAR.value) {
      otherStr = otherConst;
      for (const position of positions) {
        const posCol = table.getCol(position);
        const aarList = selectionObject[position];
        if (aarList.length === 0)
          continue;

        for (const aar of aarList) {
          aarStr = `${position} : ${aar}`;
          const splitCol = DG.Column.bool(C.COLUMNS_NAMES.SPLIT_COL, rowCount).init((i) => posCol.get(i) == aar);

          const distributionTable = DG.DataFrame.fromColumns([activityScaledCol, splitCol]);
          const currentStatsDf = model.monomerPositionStatsDf.rows.match({Pos: position, AAR: aar}).toDataFrame();
          const stats: Stats = {
            count: currentStatsDf.get(C.COLUMNS_NAMES.COUNT, 0),
            ratio: currentStatsDf.get(C.COLUMNS_NAMES.RATIO, 0),
            pValue: currentStatsDf.get(C.COLUMNS_NAMES.P_VALUE, 0),
            meanDifference: currentStatsDf.get(C.COLUMNS_NAMES.MEAN_DIFFERENCE, 0),
          };
          const distributionRoot = getDistributionAndStats(distributionTable, stats, aarStr, otherStr, true);
          $(distributionRoot).addClass('d4-flex-col');

          res.push(distributionRoot);
        }
      }
    } else if (splitByPosition.value) {
      otherStr = otherConst;
      const activityScaledData = activityScaledCol.toList();
      for (const position of positions) {
        const posCol = table.getCol(position);
        const aarList = selectionObject[position];
        if (aarList.length === 0)
          continue;

        aarStr = `${position}: {${aarList.join(', ')}}`;

        const mask = DG.BitSet.create(rowCount, (i) => aarList.includes(posCol.get(i)));
        const stats = getStats(activityScaledData, mask);
        const splitCol = DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask);
        const distributionTable = DG.DataFrame.fromColumns([activityScaledCol, splitCol]);
        const distributionRoot = getDistributionAndStats(distributionTable, stats, aarStr, otherStr, true);
        $(distributionRoot).addClass('d4-flex-col');

        res.push(distributionRoot);
      }
    } else if (splitByAAR.value) {
      const reversedSelectionObject: {[aar: string]: string[]} = {};
      const aars = [];
      for (const position of positions) {
        for (const aar of selectionObject[position]) {
          if (!reversedSelectionObject.hasOwnProperty(aar)) {
            reversedSelectionObject[aar] = [position];
            aars.push(aar);
            continue;
          }
          if (!reversedSelectionObject[aar].includes(position))
            reversedSelectionObject[aar].push(position);
        }
      }

      otherStr = otherConst;
      const activityScaledData = activityScaledCol.toList();
      for (const aar of aars) {
        const posList = reversedSelectionObject[aar];
        aarStr = `${aar}: {${posList.join(', ')}}`;

        const mask = DG.BitSet.create(rowCount, (i) => {
          const currentRow = table.row(i);
          for (const position of posList) {
            if (currentRow.get(position) == aar)
              return true;
          }
          return false;
        });
        const stats = getStats(activityScaledData, mask);
        const splitCol = DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask);
        const distributionTable = DG.DataFrame.fromColumns([activityScaledCol, splitCol]);
        const distributionRoot = getDistributionAndStats(distributionTable, stats, aarStr, otherStr, true);
        $(distributionRoot).addClass('d4-flex-col');

        res.push(distributionRoot);
      }
    } else {
      const splitCol = table.col(C.COLUMNS_NAMES.SPLIT_COL);
      if (!splitCol)
        res.push(ui.divText('No distribution'));
      else {
        otherStr = '';
        if (useSelectedStr) {
          aarStr = 'Selected';
          otherStr = otherConst;
        } else if (positionsLen) {
          aarStr = '';
          for (const position of positions) {
            const aarList = selectionObject[position];
            if (aarList.length !== 0)
              aarStr += `${position}: {${aarList.join(', ')}}; `;
          }
          if (clustersObject.length !== 0)
            aarStr += `Clusters: ${clustersObject.join(', ')}`;
          otherStr = otherConst;
        }

        const distributionTable = DG.DataFrame.fromColumns([activityScaledCol, splitCol]);
        const stats = getStats(activityScaledCol.toList(), table.selection);
        const distributionRoot = getDistributionAndStats(distributionTable, stats, aarStr, otherStr);
        $(distributionRoot).addClass('d4-flex-col');

        res.push(distributionRoot);
      }
    }
    $(distributionHost).empty().append(res);
  };

  const setDefaultProperties = (input: DG.InputBase): void => {
    input.enabled = !model.isMutationCliffSelectionEmpty;
    $(input.root).find('.ui-input-editor').css('margin', '0px');
    $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-left', '5px');
  };

  let defaultValuePos = model.splitByPos;
  let defaultValueAAR = model.splitByAAR;
  if (!model.isLogoSummarySelectionEmpty && model.isMutationCliffSelectionEmpty) {
    defaultValuePos = false;
    defaultValueAAR = false;
  }

  const splitByPosition = ui.boolInput('', defaultValuePos, updateDistributionHost);
  splitByPosition.addPostfix('Split by position');
  setDefaultProperties(splitByPosition);
  $(splitByPosition.root).css('margin-right', '10px');
  const splitByAAR = ui.boolInput('', defaultValueAAR, updateDistributionHost);
  splitByAAR.addPostfix('Split by monomer');
  setDefaultProperties(splitByAAR);

  const controlsHost = ui.divH([splitByPosition.root, splitByAAR.root]);
  const distributionHost = ui.div([], 'd4-flex-wrap');
  splitByAAR.fireChanged();

  return new DG.Widget(ui.divV([controlsHost, distributionHost]));
}

export function getDistributionAndStats(
  table: DG.DataFrame, stats: Stats, thisLabel: string, otherLabel: string = '', isTooltip: boolean = false,
): HTMLDivElement {
  const labels = ui.divV([
    ui.label(thisLabel, {style: {color: DG.Color.toHtml(otherLabel == '' ? DG.Color.blue : DG.Color.orange)}}),
    ui.label(otherLabel, {style: {color: DG.Color.toHtml(DG.Color.blue)}})]);

  const histRoot = table.plot.histogram({
    filteringEnabled: false,
    valueColumnName: table.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)?.name,
    splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
    legendVisibility: 'Never',
    showXAxis: true,
    showColumnSelector: false,
    showRangeSlider: false,
    showBinSelector: !isTooltip,
    backColor: isTooltip ? '#fdffe5' : '#fffff',
  }).root;
  histRoot.style.width = 'auto';

  const tableMap: StringDictionary = {
    'Statistics:': '',
    'Count': stats.count.toString(),
    'Ratio': stats.ratio.toFixed(2),
    'p-value': stats.pValue < 0.01 ? '<0.01' : stats.pValue.toFixed(2),
    'Mean difference': stats.meanDifference.toFixed(2),
  };

  const result = ui.divV([labels, histRoot, ui.tableFromMap(tableMap)]);
  result.style.minWidth = '200px';
  if (isTooltip)
    histRoot.style.maxHeight = '150px';
  return result;
}
