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

type Getter = [name: string, get: (look: unknown) => any];
const getters = new WeakMap<object, Getter[]>();

/** Reads the properties in one pass: every `sc.props.x` access builds wrappers for all of the viewer's properties. */
export function readLook(sc: DG.ScatterPlotViewer): ScatterLook {
  const target = sc.propertyTarget;
  const look: any = {};
  // platforms before the 2026-08 JS API have no propertyTarget
  if (target == null) {
    for (const name of LOOK_PROPERTIES)
      look[name] = (sc.props as any)[name];
    return look;
  }
  let viewerGetters = getters.get(sc.dart);
  if (!viewerGetters) {
    const byName = new Map(sc.getProperties().map((p) => [p.name, p]));
    viewerGetters = LOOK_PROPERTIES.map((name) => [name, byName.get(name)!.get]);
    getters.set(sc.dart, viewerGetters);
  }
  for (const [name, get] of viewerGetters)
    look[name] = get(target);
  return look;
}

const AXIS_TYPES = new Set<string>([DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.DATE_TIME,
  DG.COLUMN_TYPE.STRING]);

/** The core draws DOT markers with its pixel-buffer renderer unless a Markers column gives rows other shapes. */
export function isDotRenderer(look: ScatterLook): boolean {
  return look.markerType == DG.MARKER_TYPE.DOT && !look.markersColumnName;
}

/** Whether the GPU renderers draw the plot as the core would; other plots are left to the core's CPU renderers. */
export function isGpuRenderable(sc: DG.ScatterPlotViewer, look: ScatterLook): boolean {
  const df = sc.dataFrame;
  const isAxis = (name: string, map: string) => {
    const col = df.col(name);
    // a datetime axis by a time unit plots a derived column the dataframe does not hold
    return col != null && AXIS_TYPES.has(col.type) && !(col.type == DG.COLUMN_TYPE.DATE_TIME && map);
  };
  // the GPU canvas goes between the plot's canvas and its overlay; a matrix plot cell keeps the canvas only
  return sc.canvas.isConnected && sc.overlay?.parentElement === sc.canvas.parentElement &&
    look.zoomAndFilter != 'pack and zoom by filter' && !look.jitterSize && !look.jitterSizeY &&
    isAxis(look.xColumnName, look.xMap) && isAxis(look.yColumnName, look.yMap) &&
    // the canvas renderer draws labels next to the markers, the dot renderer draws none
    (isDotRenderer(look) || !look.labelColumnNames?.length || look.displayLabels == 'Never');
}
