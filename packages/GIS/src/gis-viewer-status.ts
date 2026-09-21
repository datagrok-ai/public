import * as DG from 'datagrok-api/dg';
import * as OLProj from 'ol/proj';
import {Point} from 'ol/geom';
import type {GisViewer} from './gis-viewer';

type Box = {x: number, y: number, width: number, height: number};

const MAX_POINTS = 20;

/* What the GIS map shows, for automation: its chrome and its markers as hit areas in CSS px of the
 * map canvas, and the readings a test compares. Every reading comes from the map's own view, layers
 * and feature set — never from the picture, which is base-map tiles fetched from the public
 * internet and says nothing about the product. */
export function gisViewerStatus(v: GisViewer): DG.IWidgetStatus {
  const viewport = v.root.querySelector('.ol-viewport') as HTMLElement | null;
  const canvas = viewport?.querySelector('canvas') as HTMLCanvasElement | null;
  const hitAreas: {[name: string]: Box} = {};
  const values: {[name: string]: number | string | boolean} = {};

  const origin = (canvas ?? v.root).getBoundingClientRect();
  const put = (name: string, el: Element | null | undefined) => {
    if (!el)
      return;
    const r = el.getBoundingClientRect();
    if (r.width > 0 && r.height > 0)
      hitAreas[name] = {x: r.left - origin.left, y: r.top - origin.top, width: r.width, height: r.height};
  };

  put('view', viewport);
  put('zoom in', v.root.querySelector('.ol-zoom-in'));
  put('zoom out', v.root.querySelector('.ol-zoom-out'));

  const layers = v.ol.getLayersList();
  let visible = 0;
  for (const layer of layers) {
    const name = String(layer.get('layerName') ?? '');
    if (layer.getVisible())
      visible++;
    if (name !== '')
      values[`layer "${name}" visible`] = layer.getVisible();
  }
  values['layers'] = layers.length;
  values['visible layers'] = visible;

  const panel = v.ol.panelLayersList;
  const rows = panel.dfLayersList?.rowCount ?? 0;
  if (panel.layersPanel.style.visibility !== 'hidden' && rows > 0) {
    put('layers panel', panel.layersPanel);
    const gridOrigin = panel.layersGrid.root.getBoundingClientRect();
    const dx = gridOrigin.left - origin.left;
    const dy = gridOrigin.top - origin.top;
    const cellBox = (column: string, row: number): Box | null => {
      if (!panel.layersGrid.columns.byName(column))
        return null;
      const r = panel.layersGrid.cell(column, row).bounds;
      return r.width > 0 && r.height > 0 ? {x: r.x + dx, y: r.y + dy, width: r.width, height: r.height} : null;
    };
    for (let i = 0; i < Math.min(rows, layers.length); i++) {
      const name = String(layers[i].get('layerName') ?? '');
      if (name === '')
        continue;
      const row = cellBox('name', i);
      const vis = cellBox('vis', i);
      if (row && vis)
        hitAreas[`layer "${name}"`] = {x: vis.x, y: row.y, width: row.x + row.width - vis.x, height: row.height};
      if (vis)
        hitAreas[`layer "${name}" visibility`] = vis;
    }
  }

  const view = v.ol.olMap.getView();
  const zoom = view.getZoom();
  if (zoom !== undefined)
    values['zoom'] = Math.round(zoom * 100) / 100;
  const centre = view.getCenter();
  if (centre) {
    const lonLat = OLProj.toLonLat(centre);
    values['centre'] = `${lonLat[0].toFixed(3)}, ${lonLat[1].toFixed(3)}`;
  }
  values['markers'] = v.features.length;
  values['rows shown'] = v.dataFrame.filter.trueCount;
  values['render type'] = v.renderType;

  // `getPixelFromCoordinate` answers in viewport pixels; the areas are anchored on the canvas
  const vp = viewport?.getBoundingClientRect();
  const size = Math.max(v.markerDefaultSize * 2, 4);
  let points = 0;
  for (const feature of vp ? v.features : []) {
    if (points >= MAX_POINTS)
      break;
    const geometry = feature?.getGeometry();
    if (!(geometry instanceof Point))
      continue;
    const pixel = v.ol.olMap.getPixelFromCoordinate(geometry.getCoordinates());
    if (!pixel || pixel[0] < 0 || pixel[1] < 0 || pixel[0] > vp!.width || pixel[1] > vp!.height)
      continue;
    const row = feature.get('fieldIndex');
    if (row === undefined)
      continue;
    const x = pixel[0] + vp!.left - origin.left - size / 2;
    const y = pixel[1] + vp!.top - origin.top - size / 2;
    hitAreas[`point ${row + 1}`] = {x, y, width: size, height: size};
    points++;
  }

  const parts: {[name: string]: Element} = {root: v.root};
  if (canvas)
    parts['canvas'] = canvas;
  if (hitAreas['layers panel'])
    parts['layers panel'] = v.ol.panelLayersList.layersPanel;

  return {parts, hitAreas, values, shortcuts: {}, events: [], description: null, error: null};
}
