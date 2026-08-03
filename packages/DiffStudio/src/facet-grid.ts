import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMaxGraphsInFacetGridRow} from './utils';
import {STAGE_COL_NAME} from './scripting-tools';
import {MAX_FACET_GRAPHS_COUNT} from './ui-constants';

const COLORS = DG.Color.categoricalPalette;
const COLORS_COUNT = COLORS.length;

/** Options controlling the facet grid layout. */
export interface FacetGridOptions {
  /** Independent-variable column plotted on the X axis. Defaults to the first column. */
  xColumnName?: string;
  /** Column whose categories segment-color each curve (multi-stage models). */
  segmentColumnName?: string;
  /** Upper bound on the number of plotted variables. */
  maxGraphs?: number;
}

/** Result of building a facet grid: the root element and the line-chart children it hosts. */
export interface FacetGrid {
  root: HTMLDivElement;
  plots: DG.Viewer[];
}

/** Build a faceted grid of single-variable line charts — one chart per output variable, sharing
 *  the X (argument) column. This is the single source of truth for the Diff Studio facet view,
 *  used both by the standalone app and by the `DiffStudioFacetViewer`. */
export function buildFacetGrid(table: DG.DataFrame, options?: FacetGridOptions): FacetGrid {
  const colNames = table.columns.names();

  const segmentColumnName = (options?.segmentColumnName && colNames.includes(options.segmentColumnName)) ?
    options.segmentColumnName :
    (colNames.includes(STAGE_COL_NAME) ? STAGE_COL_NAME : undefined);

  const xColumnName = (options?.xColumnName && colNames.includes(options.xColumnName)) ?
    options.xColumnName : colNames[0];

  const maxGraphs = options?.maxGraphs ?? MAX_FACET_GRAPHS_COUNT;

  const yColNames = colNames.filter((name) => name !== xColumnName && name !== segmentColumnName);
  const colsToShowCount = Math.min(maxGraphs, yColNames.length);

  if (colsToShowCount < 1)
    return {root: ui.divText('No variables to plot'), plots: []};

  const maxInRow = getMaxGraphsInFacetGridRow(colsToShowCount);

  const facetColumnPlots = new Array<DG.Viewer[]>(maxInRow);
  const plots: DG.Viewer[] = [];

  for (let i = 0; i < maxInRow; ++i)
    facetColumnPlots[i] = [];

  let idx = 0;

  for (let i = 0; i < colsToShowCount; ++i) {
    const color = COLORS[i % COLORS_COUNT];

    const plot = DG.Viewer.lineChart(table, {
      xColumnName: xColumnName,
      yColumnNames: [yColNames[i]],
      autoLayout: true,
      showXAxis: true,
      showYAxis: true,
      showXSelector: false,
      showSplitSelector: false,
      lineWidth: 2,
      segmentColumnName: segmentColumnName,
      lineColoringType: 'Custom',
      lineColor: color,
      markerColor: color,
    });

    plots.push(plot);
    facetColumnPlots[idx].push(plot);

    ++idx;

    if (idx === maxInRow)
      idx = 0;
  }

  facetColumnPlots.forEach((col) => col[col.length - 1].setOptions({showXSelector: true}));

  const rowsCount = facetColumnPlots[0].length;

  const col2Roots = (col: DG.Viewer[]) => {
    const roots = col.map((plot) => plot.root);
    if (col.length < rowsCount)
      roots.push(ui.div(''));

    return roots;
  };

  const splitCols = facetColumnPlots.map((col) => ui.splitV(col2Roots(col)));
  const root = ui.splitH(splitCols);

  return {root, plots};
} // buildFacetGrid
