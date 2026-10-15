import * as DG from 'datagrok-api/dg';

const LOOK_PROPERTIES = [
  'xColumnName', 'yColumnName', 'xAxisType', 'yAxisType', 'invertXAxis', 'invertYAxis', 'xMap', 'yMap',
  'markerType', 'markersColumnName', 'markersMap', 'sizeColumnName', 'markerOpacity', 'markerDrawBorder',
  'markerBorderWidth', 'markerMinSize', 'markerMaxSize', 'markerSizeScaling', 'showMarkersWithEmptySize',
  'colorColumnName', 'colorMap', 'colorAxisType', 'colorMin', 'colorMax', 'invertColorScheme', 'linearColorScheme',
  'categoricalColorScheme', 'filteredRowsColor', 'filteredOutRowsColor', 'selectedRowsColor', 'missingValueColor',
  'showFilteredOutPoints', 'showSelectedRows', 'rowSource',
  'zoomAndFilter', 'jitterSize', 'jitterSizeY', 'labelColumnNames', 'displayLabels',
] as const;

type LookProperty = (typeof LOOK_PROPERTIES)[number];
export type ScatterLook = Record<LookProperty, any>;

/** Marker names in the order of the core's `Marker.types`: `getMarkerTypes()` answers indexes into this list. */
export const MARKER_TYPES = [
  'circle', 'circle border', 'square', 'square border', 'cross border', 'cross x border', 'outlier',
  'diamond', 'diamond border', 'triangle top', 'triangle right', 'triangle bottom', 'triangle left',
  'triangle top border', 'triangle right border', 'triangle bottom border', 'triangle left border',
  'star', 'dot', 'gradient'];

type Getter = [name: string, get: ((look: unknown) => any) | undefined];
const getters = new WeakMap<object, Getter[]>();

/** Reads the properties in one pass: every `sc.props.x` access builds wrappers for all of the viewer's properties. */
export function readLook(sc: DG.ScatterPlotViewer): ScatterLook {
  let viewerGetters = getters.get(sc.dart);
  if (!viewerGetters) {
    const byName = new Map(sc.getProperties().map((p) => [p.name, p]));
    viewerGetters = LOOK_PROPERTIES.map((name) => [name, byName.get(name)?.get]);
    getters.set(sc.dart, viewerGetters);
  }
  // platforms before the 2026-08 JS API have no propertyTarget, and their `props` throws on a property their core
  // lacks (1.27.7: markerSizeScaling, showMarkersWithEmptySize); it reads as undefined, the behaviour before it
  const target = sc.propertyTarget;
  const look: any = {};
  for (const [name, get] of viewerGetters)
    look[name] = get && (target == null ? (sc.props as any)[name] : get(target));
  return look;
}

const AXIS_TYPES = new Set<string>([DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.DATE_TIME,
  DG.COLUMN_TYPE.STRING]);

/** The core draws DOT markers with its pixel-buffer renderer unless a Markers column gives rows other shapes. */
export function isDotRenderer(look: ScatterLook): boolean {
  return look.markerType == DG.MARKER_TYPE.DOT && !look.markersColumnName;
}

/** Why the GPU renderers cannot draw the plot as the core would, or null when they can; the plots they cannot draw
 * are left to the core's CPU renderers. */
export function whyOnCpu(sc: DG.ScatterPlotViewer, look: ScatterLook): string | null {
  const df = sc.dataFrame;
  const axisProblem = (name: string, map: string) => {
    const col = df.col(name);
    // a datetime axis by a time unit plots a derived column the dataframe does not hold
    return col == null ? `no column ${name}` : !AXIS_TYPES.has(col.type) ? `a ${col.type} axis` :
      col.type == DG.COLUMN_TYPE.DATE_TIME && map ? `a datetime axis by ${map}` : null;
  };
  // the GPU canvas goes between the plot's canvas and its overlay; a matrix plot cell keeps the canvas only
  if (!sc.canvas.isConnected || sc.overlay?.parentElement !== sc.canvas.parentElement)
    return 'not a plot of its own (a matrix plot cell, or detached)';
  if (look.zoomAndFilter == 'pack and zoom by filter')
    return 'pack and zoom by filter';
  if (look.jitterSize || look.jitterSizeY)
    return 'jitter';
  // the canvas renderer draws labels next to the markers, the dot renderer draws none
  if (!isDotRenderer(look) && look.labelColumnNames?.length && look.displayLabels != 'Never')
    return 'labels next to the markers';
  return axisProblem(look.xColumnName, look.xMap) ?? axisProblem(look.yColumnName, look.yMap);
}
