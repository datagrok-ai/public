import * as DG from 'datagrok-api/dg';
import {SARViewer, MonomerPosition, MostPotentResidues, SELECTION_MODE} from './sar-viewer';
import {COLUMNS_NAMES, TAGS} from '../utils/constants';

export function sarViewerStatus(viewer: SARViewer): DG.IWidgetStatus {
  const grid = viewer._viewerGrid;
  const visible = grid && viewer.root.contains(grid.root) && grid.canvas.getBoundingClientRect().width > 0 &&
    grid.canvas.getBoundingClientRect().height > 0;
  const inner = visible ? grid.getWidgetStatus() : null;
  const hitAreas: DG.IWidgetStatus['hitAreas'] = {};
  const values: NonNullable<DG.IWidgetStatus['values']> = {
    'activity scaling': viewer.activityScaling,
    'data source': viewer.dataSource,
  };
  const status: DG.IWidgetStatus = {
    parts: {root: viewer.root, ...inner?.parts},
    hitAreas, values, shortcuts: {}, events: [], description: null,
    error: viewer.root.textContent?.includes('Please, select a sequence and activity columns') ?
      'Please, select a sequence and activity columns in the viewer properties' : null,
  };
  if (viewer.dataFrame)
    values['rows shown'] = viewer.dataSource === 'All' ? viewer.dataFrame.rowCount : viewer.dataFrame.filter.trueCount;
  if (viewer._positionColumns)
    values.positions = viewer._positionColumns.length;

  const stats = viewer._monomerPositionStats;
  if (stats) {
    let cells = 0;
    const monomers = new Set<string>();
    for (const [position, positionStats] of Object.entries(stats)) {
      if (position === 'general')
        continue;
      for (const [monomer, cell] of Object.entries(positionStats)) {
        if (monomer === 'general' || !('count' in cell) || cell.count <= 0)
          continue;
        cells++;
        monomers.add(monomer);
      }
    }
    values['cells with stats'] = cells;
    values.monomers = monomers.size;
  }

  const cliffs = viewer._mutationCliffs;
  if (cliffs) {
    let cells = 0;
    let pairs = 0;
    const uniquePairs = new Set<string>();
    for (const positions of cliffs.values()) {
      for (const indexes of positions.values()) {
        if (indexes.size > 0)
          cells++;
        for (const [from, targets] of indexes) {
          pairs += targets.length;
          for (const to of targets)
            uniquePairs.add(from < to ? `${from}-${to}` : `${to}-${from}`);
        }
      }
    }
    values['cliff cells'] = cells;
    values['cliff pairs'] = pairs;
    values['unique cliff pairs'] = uniquePairs.size;
  }

  const mode = viewer instanceof MonomerPosition && viewer.dataFrame ? viewer.mode : SELECTION_MODE.INVARIANT_MAP;
  if (viewer instanceof MonomerPosition)
    values.mode = mode;
  const selection = mode === SELECTION_MODE.MUTATION_CLIFFS ?
    viewer._mutationCliffsSelection : viewer._invariantMapSelection;
  values['selected monomer-positions'] = Object.entries(selection ?? {})
    .flatMap(([position, monomers]) => monomers.map((monomer) => `${position}:${monomer}`)).join(', ');

  if (!inner || !grid)
    return status;

  const canvasBounds = grid.canvas.getBoundingClientRect();
  hitAreas.view = {x: 0, y: 0, width: canvasBounds.width, height: canvasBounds.height};
  status.shortcuts.ContextMenu = 'view';
  for (const [name, bounds] of Object.entries(inner.hitAreas)) {
    if (name.startsWith('header ') || name.includes('scroll '))
      hitAreas[name] = bounds;
    const cell = /^cell (\d+) of (.+)$/.exec(name);
    if (!cell)
      continue;
    const row = Number(cell[1]) - 1;
    const column = cell[2];
    const monomer = grid.dataFrame.get(COLUMNS_NAMES.MONOMER, row) as string;
    if (viewer instanceof MonomerPosition) {
      if (column === COLUMNS_NAMES.MONOMER)
        hitAreas[`monomer ${monomer}`] = bounds;
      const cellStats = stats?.[column]?.[monomer];
      if (!viewer._positionColumns?.some((position) => position.name === column) || !cellStats?.count)
        continue;
      const key = `cell ${monomer} at ${column}`;
      hitAreas[key] = bounds;
      values[`count of ${key}`] = cellStats.count;
      values[`mean difference of ${key}`] = cellStats.meanDifference;
      values[`p-value of ${key}`] = cellStats.pValue ?? '';
      values[`value of ${key}`] = cellStats.aggValue ?? cellStats.count;
      if (cliffs) {
        const indexes = cliffs.get(monomer)?.get(column);
        values[`cliffs of ${key}`] = indexes ? Array.from(indexes.values())
          .reduce((sum, targets) => sum + targets.length, 0) : 0;
      }
      if (mode === SELECTION_MODE.INVARIANT_MAP) {
        const position = viewer._positionColumns.find((position) => position.name === column)!;
        const color = position.temp[TAGS.INVARIANT_MAP_COLOR_CACHE]?.[monomer];
        if (color !== undefined)
          values[`color of ${key}`] = DG.Color.toHtml(color);
      }
    } else if (viewer instanceof MostPotentResidues) {
      const position = grid.dataFrame.get(COLUMNS_NAMES.POSITION, row);
      if (column === 'Diff')
        hitAreas[`position ${position}`] = bounds;
      else if (column === COLUMNS_NAMES.MONOMER)
        hitAreas[`monomer of position ${position}`] = bounds;
    }
  }

  if (viewer instanceof MonomerPosition && viewer.monomerSearchInput.input.getBoundingClientRect().width > 0) {
    const bounds = viewer.monomerSearchInput.input.getBoundingClientRect();
    hitAreas.search = {x: bounds.left - canvasBounds.left, y: bounds.top - canvasBounds.top,
      width: bounds.width, height: bounds.height};
  }
  if (viewer instanceof MostPotentResidues) {
    const table = grid.dataFrame;
    values.positions = table.rowCount;
    values['rows shown'] = table.filter.trueCount;
    const columns = {
      'monomer': COLUMNS_NAMES.MONOMER,
      'mean difference': COLUMNS_NAMES.MEAN_DIFFERENCE,
      'p-value': COLUMNS_NAMES.P_VALUE,
      'count': COLUMNS_NAMES.COUNT,
      'ratio': COLUMNS_NAMES.RATIO,
    };
    for (let row = 0; row < table.rowCount; row++) {
      const position = table.get(COLUMNS_NAMES.POSITION, row);
      for (const [reading, name] of Object.entries(columns)) {
        const column = table.getCol(name);
        values[`${reading} at ${position}`] = column.isNone(row) ? '' : column.get(row);
      }
    }
  }
  return status;
}
