import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pmpo.css';

import {COLORS, DESCR_TABLE_TITLE, DESCR_TITLE, DescriptorStatistics, P_VAL, SELECTED,
  STAT_TO_TITLE_MAP} from './pmpo-defs';

export function getDescriptorStatisticsTable(stats: Map<string, DescriptorStatistics>): DG.DataFrame {
  const descrCount = stats.size;
  const rawArrs = new Map<string, Float64Array>();

  // Create raw data arrays
  STAT_TO_TITLE_MAP.forEach((_, key) => {
    rawArrs.set(key, new Float64Array(descrCount));
  });

  const descrNames = [...stats.keys()];
  const cols = [DG.Column.fromStrings(DESCR_TITLE, descrNames)];

  // Fill stat columns
  descrNames.forEach((descr, idx) => {
    const curStat = stats.get(descr);

    if (curStat != null) {
      STAT_TO_TITLE_MAP.forEach((_, key) => {
        const val = curStat[key as keyof DescriptorStatistics];
        const arr = rawArrs.get(key);
        arr![idx] = val;
      });
    }
  });

  // Create stat columns
  STAT_TO_TITLE_MAP.forEach((title, field) => {
    cols.push(DG.Column.fromFloat64Array(title, rawArrs.get(field)!));
  });

  // Create the resulting table
  const res = DG.DataFrame.fromColumns(cols);
  res.name = DESCR_TABLE_TITLE;

  return res;
} // getDescriptorStatisticsTable

export function getFilteredByPvalue(descrStats: DG.DataFrame, pValThresh: number): string[] {
  const selected: string[] = [];

  const descrCol = descrStats.col(DESCR_TITLE);

  if (descrCol == null)
    throw new Error(`No column "${DESCR_TITLE} in the table with descriptors statistics.`);

  const descr = descrCol.toList();

  const pValCol = descrStats.col(P_VAL);

  if (pValCol == null)
    throw new Error(`No column "${P_VAL} in the table with descriptors statistics.`);

  const pVals = pValCol.getRawData();

  for (let i = 0; i < descrStats.rowCount; ++i) {
    if (pVals[i] < pValThresh)
      selected.push(descr[i]);
  }

  return selected;
} // getFilteredByPvalue

export function addSelectedDescriptorsCol(descrStats: DG.DataFrame, selected: string[]): DG.DataFrame {
  if (selected.length < 1)
    throw new Error('Empty list of selected descriptors.');

  const rowCount = descrStats.rowCount;
  const selArr = new Array<boolean>(rowCount);
  const descrCol = descrStats.col(DESCR_TITLE);

  if (descrCol == null)
    throw new Error(`No column "${DESCR_TITLE} in the table with descriptors statistics.`);

  const descr = descrCol.toList();
  let res = true;
  const colors: Record<string, string> = {};

  for (let i = 0; i < rowCount; ++i) {
    res = selected.includes(descr[i]);
    selArr[i] = res;
    colors[descr[i]] = res ? COLORS.SELECTED : COLORS.SKIPPED;
  }

  descrCol.colors.setCategorical(colors);

  // Added selected column
  descrStats.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, SELECTED, selArr));

  return descrStats;
} // addSelectedDescriptorsCol

export function getDescriptorStatisticsGrid(table: DG.DataFrame,
  selectedByPvalue: string[], selectedByCorr: string[]): DG.Grid {
  const grid = DG.Viewer.grid(table, {
    showTitle: true,
    title: table.name,
  });

  grid.sort([SELECTED, P_VAL], [false, true]);
  grid.col(P_VAL)!.format = 'scientific';

  // set tooltips
  grid.onCellTooltip(function(cell, x, y) {
    if (cell.isColHeader) {
      const cellCol = cell.tableColumn;
      if (cellCol) {
        if (cell.tableColumn.name === DESCR_TITLE) {
          ui.tooltip.show(getDescrTooltip(), x, y);

          return true;
        }

        return false;
      }
    } else {
      if (cell.isTableCell) {
        const cellCol = cell.tableColumn;
        if (cellCol) {
          if (cell.tableColumn.name === DESCR_TITLE) {
            const value = cell.value;
            if (selectedByCorr.includes(value))
              ui.tooltip.show('Selected for model construction.', x, y);
            else if (selectedByPvalue.includes(value))
              ui.tooltip.show('Excluded due to a high correlation with other descriptors.', x, y);
            else
              ui.tooltip.show('Excluded due to a high p-value.', x, y);

            return true;
          }

          return false;
        }
      }
    }
  });

  return grid;
} // getDescriptorStatisticsGrid

function getDescrTooltip(): HTMLElement {
  const firstLine = ui.div();
  firstLine.classList.add('eda-pmpo-tooltip-line');
  const selectedBox = ui.div();
  selectedBox.classList.add('eda-pmpo-box');
  selectedBox.style.backgroundColor = COLORS.SELECTED;
  const selectedLabel = ui.span([]);
  selectedLabel.textContent = '- selected';
  firstLine.appendChild(selectedBox);
  firstLine.appendChild(selectedLabel);

  const secondLine = ui.div();
  secondLine.classList.add('eda-pmpo-tooltip-line');
  const nonSelectedBox = ui.div();
  nonSelectedBox.classList.add('eda-pmpo-box');
  nonSelectedBox.style.backgroundColor = COLORS.SKIPPED;
  const nonSelectedLabel = ui.span([]);
  nonSelectedLabel.textContent = '- excluded';

  secondLine.appendChild(nonSelectedBox);
  secondLine.appendChild(nonSelectedLabel);

  return ui.divV([
    ui.h2(DESCR_TITLE),
    ui.divText('Use of descriptors in model construction:'),
    firstLine,
    secondLine,
  ]);
} // getDescrTooltip

export function getFilteredByCorrelations(descriptors: DG.ColumnList, selectedByPvalue: string[],
  statistics: Map<string, DescriptorStatistics>, r2Tresh: number): string[] {
  const getCorrelations = () => {
    const triples: [string, string, number][] = [];

    for (let i = 0; i < selectedByPvalue.length; ++i) {
      for (let j = i + 1; j < selectedByPvalue.length; ++j) {
        triples.push([
          selectedByPvalue[i],
          selectedByPvalue[j],
          descriptors.byName(selectedByPvalue[i]).stats.corr(descriptors.byName(selectedByPvalue[j]))**2,
        ]);
      }
    }

    return triples;
  };

  const correlations = getCorrelations().sort((a, b) => b[2] - a[2]);

  const keep = new Set(selectedByPvalue);

  correlations.filter((triple) => triple[2] > r2Tresh).forEach((triple) => {
    const [descr1, descr2, _] = triple;
    const pVal1 = statistics.get(descr1)!.pValue;
    const pVal2 = statistics.get(descr2)!.pValue;
    const tStat1 = statistics.get(descr1)!.tstat;
    const tStat2 = statistics.get(descr2)!.tstat;

    if (pVal1 > pVal2)
      keep.delete(descr1);
    else if (pVal1 < pVal2)
      keep.delete(descr2);
    else { // process the case of p-value = 0 (or too small)
      if (Math.abs(tStat1) > Math.abs(tStat2))
        keep.delete(descr2);
      else
        keep.delete(descr1);
    }
  });

  console.log(keep);

  return [...keep];
} // getFilteredByCorrelations
