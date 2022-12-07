/* eslint-disable block-spacing */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {OpenLayers} from '../src/gis-openlayer';

export const SEMTYPEGIS = {
  LONGITUDE: 'Longitude',
  LATIITUDE: 'Latitude',
  ALTITUDE: 'gis-altitude',
  GISCOUNTRY: 'gis-country',
  GISSTATE: 'gis-state',
  GISADDRESS: 'gis-address',
  GISZIPCODE: 'gis-zipcode',

  GISOBJECT: 'gis-object',
  GISPOINT: 'gis-point',
  GISAREA: 'gis-area',
};

export type gisCoordinate = [number, number, number?];
export type gisFeatureProperties = {[x: string]: any};
export type gisPolygonCoords = Array<gisCoordinate>;
export type gisPolygon = Array<gisPolygonCoords>;
export type gisPolygons = Array<gisPolygon>;

export class GisArea {
  semtype: string = SEMTYPEGIS.GISAREA;
  // coordinates: Array<gisCoordinate>;
  coordinates: gisPolygons;
  attributes: gisFeatureProperties = {};
  // maprefference: OpenLayers | null; //TODO: maybe we should refference to superclass (e.g. MapEngine above OpenLayers)

  constructor(coord: gisPolygons, attr?: gisFeatureProperties, parentmap: OpenLayers | null = null) {
    this.coordinates = coord;
    if (attr)
      this.attributes = attr;
    // this.maprefference = parentmap;
  }

  toString() {
    let strRes = '';
    strRes = JSON.stringify(this);
    return strRes;
  }
}

function drawContourByCoords(g: CanvasRenderingContext2D,
  x: number, y: number, w: number, h: number,
  polygonsarray: gisPolygons) {
  //coordinates for area polygons are stored in arrays like this
  //[[polygon1:[outer contour:[...]], [inner contour:[...]],...], polygon1:[outer contour:[...]]]
  if (polygonsarray.length === 0)
    return;
  //detect scale>>
  let xMin = Number.MAX_SAFE_INTEGER;
  let xMax = Number.MIN_SAFE_INTEGER;
  let yMin = Number.MAX_SAFE_INTEGER;
  let yMax = Number.MIN_SAFE_INTEGER;
  let xScale = 1;
  let yScale = 1;
  for (let p = 0; p < polygonsarray.length; p++) {
    const parray: gisPolygon = polygonsarray[p];
    for (let i = 0; i < parray.length; i++) {
      const coordarray: gisPolygonCoords = parray[i];
      for (let k = 0; k < coordarray.length - 1; k++) {
        for (let i = 0; i < coordarray.length; i++) {
          xMin = Math.min(xMin, coordarray[i][0]);
          xMax = Math.max(xMax, coordarray[i][0]);
          yMin = Math.min(yMin, coordarray[i][1]);
          yMax = Math.max(yMax, coordarray[i][1]);
        }
        xScale = (w - 10) / Math.abs(xMax - xMin);
        yScale = (h - 10) / Math.abs(yMax - yMin);
        if (yScale < xScale)
          xScale = yScale;
        else yScale = xScale;
      }
    }
  }

  //centering coefficients of contour for canvas
  const xShift = ((w - 10) - Math.abs((xMax - xMin) * xScale)) / 2;
  const yShift = ((h - 10) - Math.abs((yMax - yMin) * yScale)) / 2;
  //draw contour>>
  for (let p = 0; p < polygonsarray.length; p++) {
    const parray: gisPolygon = polygonsarray[p];
    for (let i = 0; i < parray.length; i++) {
      g.fillStyle = '#FEEEEE';
      g.strokeStyle = '#0000F0';
      g.beginPath();
      const coordarray: gisPolygonCoords = parray[i];
      for (let k = 0; k < coordarray.length - 1; k++) {
        const x1 = (x + 5 + xShift) + (coordarray[k][0] - xMin) * xScale;
        const y1 = (y + h - 5 - yShift) - (coordarray[k][1] - yMin) * yScale;
        const x2 = (x + 5 + xShift) + (coordarray[k+1][0] - xMin) * xScale;
        const y2 = (y + h - 5 - yShift) - (coordarray[k+1][1] - yMin) * yScale;
        g.moveTo(x1, y1);
        g.lineTo(x2, y2);
      } //<<for contour coordinates
      g.closePath();
      g.stroke();
    } //<<for polygon contours
  } //<<for polygons
}

export class GisAreaCanvasRenderer extends DG.CanvasRenderer {
  get defaultWidth(): number | null { return 200; }
  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    obj: GisArea, context: any): void {
    drawContourByCoords(g, x, y, w, h, obj.coordinates);
  }
}

//name: gisAreaGridCellRenderer
//tags: cellRenderer
//meta.cellType: gis-area
//output: grid_cell_renderer result
export function gisAreaGridCellRenderer() {
  return new GisAreaGridCellRenderer();
}

export class GisAreaGridCellRenderer extends DG.GridCellRenderer {
  constructor() {
    super();
  }
  get cellType() { return SEMTYPEGIS.GISAREA; }
  get defaultWidth() { return 100; }
  get defaultHeight() { return 50; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ): void {
    const cellVal = gridCell.cell.value;
    let cellObj: any;
    if ((cellVal instanceof GisArea))
      cellObj = cellVal;
    if (typeof cellVal === 'string')
      cellObj = JSON.parse(cellVal);

    // const objArea = new GisArea(cellObj.coordinates, cellObj.attributes);
    if (cellObj.coordinates)
      drawContourByCoords(g, x, y, w, h, cellObj.coordinates);
  }
  //<< end of GisAreaGridCellRenderer class
}

//Area semantic type handler
export class GisAreaHandler extends DG.ObjectHandler {
  get type() { return SEMTYPEGIS.GISAREA; }
  get name() { return SEMTYPEGIS.GISAREA + ' handler'; }

  isApplicable(obj: any) {
    // console.log('isApplicable ' + obj + ' = ' + (obj instanceof GisArea));
    if (obj.hasOwnProperty('semType')) { //TODO: check why we cant find semType in SemanticValue object
      console.log(obj + ' semType = ' + (obj.semType));
      if (obj.semType === SEMTYPEGIS.GISAREA) return true;
    }
    if (obj.hasOwnProperty('value')) { //TODO: check why we cant find value in SemanticValue object
      console.log(obj + ' value = ' + (obj.value));
      if (obj.value instanceof GisArea) return true;
    }

    return (obj instanceof GisArea);
  }

  getCanvasRenderer() { return new GisAreaCanvasRenderer(); }
  getGridCellRenderer() { return new GisAreaGridCellRenderer(); }

  renderIcon() { return ui.iconFA('bullseye'); }
  renderTooltip(obj: GisArea) { return ui.divText(`Area of [${obj.coordinates.length}] vertexes`); }

  renderProperties(obj: GisArea) {
    //prepare subelements
    const coordsDF = DG.DataFrame.fromObjects(obj.coordinates);
    let grid: DG.Viewer;
    if (coordsDF)
      grid = DG.Viewer.fromType('Grid', coordsDF);

    const acc = ui.accordion();
    acc.addPane('Details', () => ui.tableFromMap(obj.attributes));
    acc.addPane('Coordinates', () => {
      if (!grid)
        return ui.div('No coordinates');
      grid.root.style.maxWidth = '300px';
      return ui.div(grid.root);
    });
    acc.addPane('Contour', () => {
      const canvas = ui.canvas(200, 100);
      this.getCanvasRenderer().render(canvas.getContext('2d')!, 0, 0, 200, 100, obj, null);
      return ui.div([canvas]);
    });

    const panelProperties = ui.divV([ui.divText(`Properties for GisObject`), acc.root]);
    return panelProperties;
  }

  renderMarkup(obj: GisArea) {
    const m = ui.span([this.renderIcon(), ui.label('Area')]);
    m.style.color = 'red';
    return m;
  }

  renderCard(obj: GisArea, context: any) {
    return ui.bind(obj, ui.divV([
      this.renderMarkup(obj),
      ui.divText(`Context: ${context}`)], 'd4-gallery-item'));
  }
}

/*
//GisPoint class
*/
export class GisPoint {
  semtype: string = SEMTYPEGIS.GISPOINT;
  coordinates: gisCoordinate = [0, 0, 0];
  attributes: gisFeatureProperties = {};
  // attributes: Record<string, any> = {};

  // constructor(lng: number, lat: number, alt?: number, attr?: Record<string, any>) {
  constructor(lng: number, lat: number, alt?: number, attr?: gisFeatureProperties) {
    this.coordinates[0] = lng;
    this.coordinates[1] = lat;
    if (alt)
      this.coordinates[2] = alt;
    if (attr)
      this.attributes = attr;
  }
  get x(): number { return this.coordinates[0]; }
  get y(): number { return this.coordinates[1]; }
  get z(): number { return this.coordinates[2] ? this.coordinates[2] : 0; }
  set x(val: number) { this.coordinates[0] = val; }
  set y(val: number) { this.coordinates[1] = val; }
  set z(val: number) { this.coordinates[2] = val; }
  get lng(): number { return this.coordinates[0]; }
  get lat(): number { return this.coordinates[1]; }
  get alt(): number { return this.coordinates[2] ? this.coordinates[2] : 0; }
  set lng(val: number) { this.coordinates[0] = val; }
  set lat(val: number) { this.coordinates[1] = val; }
  set alt(val: number) { this.coordinates[2] = val; }

  toString() {
    let strRes = '';
    strRes = JSON.stringify(this);
    return strRes;
    // return 'Gis Point Object';
  }
}

export class GisPointGridCellRenderer extends DG.GridCellRenderer {
  constructor() {
    super();
  }

  get cellType() { return SEMTYPEGIS.GISPOINT; }
  get defaultWidth() { return 100; }
  get defaultHeight() { return 50; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ): void {
    g.fillStyle = 'darkgray';
    gridCell.customText = 'GisPoint!';
    // gridCell.ob
    // g.fillText('['+this.longitude+';'+this.latitude+']', x + 5, y + 5);
    g.fillText('[GISPoint]', x + 5, y + 5);
  }
}

export class GisPointCanvasRenderer extends DG.CanvasRenderer {
  get defaultWidth(): number | null {
    return null;
  }
  get defaultHeight(): number | null {
    return null;
  }
  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    obj: GisPoint, context: any) {
    g.fillText('[' + obj.x + ';' + obj.y + ']', x + 10, y + 10);
  }
}

//Point semantic type handler
export class GisPointHandler extends DG.ObjectHandler {
  get type() { return SEMTYPEGIS.GISPOINT; }

  isApplicable(obj: any) { return (obj instanceof GisPoint); }

  // getCanvasRenderer() {return null;}
  getCanvasRenderer() { return new GisPointCanvasRenderer(); }
  getGridCellRenderer() { return new GisPointGridCellRenderer(); }

  renderIcon() { return ui.iconFA('bullseye'); }
  renderTooltip(obj: GisPoint) { return ui.divText(`[${obj.x} ; ${obj.y}]`); }
  renderProperties(obj: GisPoint) {
    const panelProperties = ui.divV([]);
    panelProperties.appendChild(ui.divText(`Properties for Point`));
    panelProperties.appendChild(ui.divText(`[${obj.x} ; ${obj.y}]`));
    return panelProperties;
  }

  renderMarkup(obj: GisPoint) {
    const m = ui.span([this.renderIcon(), ui.label('Point')]);
    m.style.color = 'red';
    return m;
  }

  renderCard(obj: GisPoint, context: any) {
    return ui.bind(obj, ui.divV([
      this.renderMarkup(obj),
      ui.divText(`Context: ${context}`)], 'd4-gallery-item'));
  }
}
