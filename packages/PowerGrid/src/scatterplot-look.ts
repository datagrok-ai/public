import * as DG from 'datagrok-api/dg';

const LOOK_PROPERTIES = [
  'xColumnName', 'yColumnName', 'xAxisType', 'yAxisType', 'invertXAxis', 'invertYAxis', 'xMap', 'yMap',
  'markerType', 'markersColumnName', 'markersMap', 'sizeColumnName', 'markerOpacity', 'markerDrawBorder',
  'markerBorderWidth', 'markerMinSize', 'markerMaxSize', 'markerSizeScaling', 'showMarkersWithEmptySize',
  'colorColumnName', 'colorMap', 'colorAxisType', 'colorMin', 'colorMax', 'invertColorScheme', 'linearColorScheme',
  'categoricalColorScheme', 'filteredRowsColor', 'filteredOutRowsColor', 'selectedRowsColor', 'missingValueColor',
  'showFilteredOutPoints', 'showSelectedRows', 'rowSource',
  'zoomAndFilter', 'jitterSize', 'jitterSizeY', 'labelColumnNames', 'displayLabels', 'xCategoryOrder', 'yCategoryOrder',
] as const;

/** What deciding between the GPU and the CPU reads: the core asks for it twice on every mouse move. */
const CHECK_PROPERTIES = [
  'xColumnName', 'yColumnName', 'xMap', 'yMap', 'markerType', 'markersColumnName', 'sizeColumnName', 'markerMaxSize',
  'markerBorderWidth', 'zoomAndFilter', 'jitterSize', 'jitterSizeY', 'labelColumnNames', 'displayLabels',
  'xCategoryOrder', 'yCategoryOrder',
] as const satisfies readonly LookProperty[];

type LookProperty = (typeof LOOK_PROPERTIES)[number];
export type ScatterLook = Record<LookProperty, any>;
export type CheckLook = Pick<ScatterLook, (typeof CHECK_PROPERTIES)[number]>;

/** Marker names in the order of the core's `Marker.types`: `getMarkerTypes()` answers indexes into this list. */
export const MARKER_TYPES = [
  'circle', 'circle border', 'square', 'square border', 'cross border', 'cross x border', 'outlier',
  'diamond', 'diamond border', 'triangle top', 'triangle right', 'triangle bottom', 'triangle left',
  'triangle top border', 'triangle right border', 'triangle bottom border', 'triangle left border',
  'star', 'dot', 'gradient'];

type Getter = (look: unknown) => any;
const getters = new WeakMap<object, Map<string, Getter | undefined>>();

/** Reads the properties in one pass: every `sc.props.x` access builds wrappers for all of the viewer's properties. */
function read<K extends LookProperty>(sc: DG.ScatterPlotViewer, names: readonly K[]): Pick<ScatterLook, K> {
  let viewerGetters = getters.get(sc.dart);
  if (!viewerGetters) {
    const byName = new Map(sc.getProperties().map((p) => [p.name, p]));
    viewerGetters = new Map(LOOK_PROPERTIES.map((name) => [name, byName.get(name)?.get]));
    getters.set(sc.dart, viewerGetters);
  }
  // platforms before the 2026-08 JS API have no propertyTarget, and their `props` throws on a property their core
  // lacks (1.27.7: markerSizeScaling, showMarkersWithEmptySize); it reads as undefined, the behaviour before it
  const target = sc.propertyTarget;
  const look: any = {};
  for (const name of names) {
    const get = viewerGetters.get(name);
    look[name] = get && (target == null ? (sc.props as any)[name] : get(target));
  }
  return look;
}

export function readLook(sc: DG.ScatterPlotViewer): ScatterLook {
  return read(sc, LOOK_PROPERTIES);
}

export function readCheckLook(sc: DG.ScatterPlotViewer): CheckLook {
  return read(sc, CHECK_PROPERTIES);
}

const canvases = new WeakMap<object, {canvas: HTMLCanvasElement, overlay: HTMLCanvasElement}>();

/** The plot's canvas and overlay. The core creates both once per viewer, while `sc.canvas` and `sc.overlay` build
 * the viewer's whole info map, column selectors and legend included, on every access. */
export function plotCanvases(sc: DG.ScatterPlotViewer): {canvas: HTMLCanvasElement, overlay: HTMLCanvasElement} {
  let elements = canvases.get(sc.dart);
  if (!elements) {
    const info = sc.getInfo();
    canvases.set(sc.dart, elements = {canvas: info['canvas'], overlay: info['overlay']});
  }
  return elements;
}

const AXIS_TYPES = new Set<string>([DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.DATE_TIME,
  DG.COLUMN_TYPE.STRING]);

/** The core draws DOT markers with its pixel-buffer renderer unless a Markers column gives rows other shapes. */
export function isDotRenderer(look: Pick<ScatterLook, 'markerType' | 'markersColumnName'>): boolean {
  return look.markerType == DG.MARKER_TYPE.DOT && !look.markersColumnName;
}

/** Whether the GPU renderers draw the plot as the core would; other plots are left to the core's CPU renderers. */
export function isGpuRenderable(sc: DG.ScatterPlotViewer, look: CheckLook): boolean {
  const df = sc.dataFrame;
  const isAxis = (name: string, map: string) => {
    const col = df.col(name);
    // a datetime axis by a time unit plots a derived column the dataframe does not hold
    return col != null && AXIS_TYPES.has(col.type) && !(col.type == DG.COLUMN_TYPE.DATE_TIME && map);
  };
  const isOrdered = (name: string, order: string) => !!order && df.col(name)?.type == DG.COLUMN_TYPE.STRING;
  // the GPU canvas goes between the plot's canvas and its overlay; a matrix plot cell keeps the canvas only
  const {canvas, overlay} = plotCanvases(sc);
  return canvas.isConnected && overlay?.parentElement === canvas.parentElement &&
    look.zoomAndFilter != 'pack and zoom by filter' && !look.jitterSize && !look.jitterSizeY &&
    // a category order moves the categories to slots the GPU data does not hold; it is kept, unused, on a numeric axis
    !isOrdered(look.xColumnName, look.xCategoryOrder) && !isOrdered(look.yColumnName, look.yCategoryOrder) &&
    isAxis(look.xColumnName, look.xMap) && isAxis(look.yColumnName, look.yMap) &&
    // the canvas renderer draws labels next to the markers, the dot renderer draws none
    (isDotRenderer(look) || !look.labelColumnNames?.length || look.displayLabels == 'Never');
}
