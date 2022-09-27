/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const SEMTYPEGIS = {
  LONGITUDE: 'gis-longitude',
  LATIITUDE: 'gis-latitude',
  ALTITUDE: 'gis-altitude',
  GISPOINT: 'gis-point',
  GISAREA: 'gis-area',
  GISOBJECT: 'gis-object',
  GISCOUNTRY: 'gis-country',
  GISSTATE: 'gis-state',
  GISADDRESS: 'gis-address',
  GISZIPCODE: 'gis-zipcode',
};

export type gisCoordinate = [number, number, number?];
export type gisFeatureProperties = {[x: string]: any};

export class GisArea {
  semtype: string = SEMTYPEGIS.GISAREA;
  coordinates: Array<gisCoordinate>; //Array<[number, number, number?]>
  attributes: gisFeatureProperties = {};
  constructor(coord: Array<gisCoordinate>, attr?: gisFeatureProperties) {
    this.coordinates = coord;
    if (attr) this.attributes = attr;
  }

  toString() {
    let strRes = '';
    strRes = JSON.stringify(this);
    return strRes;
    // return 'Gis Area Object';
  }
}

export class GisAreaCanvasRenderer extends DG.CanvasRenderer {
  get defaultWidth(): number | null {
    return 200; //null;
  }
  get defaultHeight(): number | null {
    return 100; //null;
  }
  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    obj: GisArea, context: any): void {
    if (obj.coordinates.length == 0) return;
    //detect scale>>
    let xMin = obj.coordinates[0][0];
    let xMax = obj.coordinates[0][0];
    let yMin = obj.coordinates[0][1];
    let yMax = obj.coordinates[0][1];
    for (let i = 0; i < obj.coordinates.length; i++) {
      if (obj.coordinates[i][0] < xMin) xMin = obj.coordinates[i][0];
      if (obj.coordinates[i][0] > xMax) xMax = obj.coordinates[i][0];
      if (obj.coordinates[i][1] < yMin) yMin = obj.coordinates[i][1];
      if (obj.coordinates[i][1] > yMax) yMax = obj.coordinates[i][1];
    }
    let xScale = (w-10)/Math.abs(xMax-xMin);
    let yScale = (h-10)/Math.abs(yMax-yMin);
    if (yScale < xScale) xScale = yScale;
    else yScale = xScale;
    //draw contour>>
    g.fillStyle = '#FEEEEE';
    g.strokeStyle = '#FF0000';
    g.beginPath();
    for (let i = 0; i < obj.coordinates.length-1; i++) {
      const x1 = (x+5)+(obj.coordinates[i][0]-xMin)*xScale;
      const y1 = (h-5)-(obj.coordinates[i][1]-yMin)*yScale;
      const x2 = (x+5)+(obj.coordinates[i+1][0]-xMin)*xScale;
      const y2 = (h-5)-(obj.coordinates[i+1][1]-yMin)*yScale;
      g.moveTo(x1, y1);
      g.lineTo(x2, y2);
    }
    g.closePath();
    g.stroke();
  }
}

//name: gisAreaWidget
//tags: panel, widgets
//input: object gisArea {semType: gis-area}
//output: widget result
//condition: true
export function gisAreaWidget(gisArea: any): DG.Widget | null {
//this is temporary code - should be filled with usefull functionality
  if ((!gisArea) && !(gisArea instanceof GisArea)) return null;

  const strToAdd: string = (gisArea as GisArea).semtype;
  let widgetStyle: DG.ElementOptions = { };
  widgetStyle = {style: {'color': '#F55'}};

  return new DG.Widget(ui.divText('gis Area widget ' + strToAdd, widgetStyle));
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

  get cellType() {return SEMTYPEGIS.GISAREA;}
  get defaultWidth() {return 100;}
  get defaultHeight() {return 50;}

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ): void {
    g.fillStyle = 'darkgray';
    gridCell.customText = 'GisArea!';
    // gridCell.ob
    // g.fillText('['+this.longitude+';'+this.latitude+']', x + 5, y + 5);
    g.fillText('[GISArea]', x + 5, y - 5);
  }
}

//Area semantic type handler
export class GisAreaHandler extends DG.ObjectHandler {
  get type() {return SEMTYPEGIS.GISAREA;}
  get name() {return SEMTYPEGIS.GISAREA+ ' handler';}

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

  // getCanvasRenderer() {return null;}
  getCanvasRenderer() {return new GisAreaCanvasRenderer();}
  getGridCellRenderer() {return new GisAreaGridCellRenderer();}

  renderIcon() {return ui.iconFA('bullseye');}
  renderTooltip(obj: GisArea) {return ui.divText(`Area of [${obj.coordinates.length}] vertexes`);}

  renderProperties(obj: GisArea) {
    //prepare subelements
    const coordsDF = DG.DataFrame.fromObjects(obj.coordinates);
    let grid: DG.Viewer;
    if (coordsDF) grid = DG.Viewer.fromType('Grid', coordsDF);

    const acc = ui.accordion();
    acc.addPane('Details', () => ui.tableFromMap(obj.attributes));
    acc.addPane('Coordinates', () => {
      if (!grid) return ui.div('No coordinates');
      grid.root.style.maxWidth = '300px';
      return ui.div(grid.root);
    });
    acc.addPane('Contour', () => {
      const canvas = ui.canvas(200, 100);
      this.getCanvasRenderer().render(canvas.getContext('2d')!, 0, 0, 200, 100, obj, null);
      return ui.div([canvas]);
    });

    const panelProperties = ui.divV([ui.divText(`Properties for Area`), acc.root]);
    // panelProperties.appendChild(ui.divText(`Properties for Area`));
    // for (let i = 0; i < obj.coordinates.length; i++)
    //   panelProperties.appendChild(ui.divText(`[${obj.coordinates[i][0]} ; ${obj.coordinates[i][1]}]`));
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
    if (alt) this.coordinates[2] = alt;
    if (attr) this.attributes = attr;
  }
  get x(): number {return this.coordinates[0];}
  get y(): number {return this.coordinates[1];}
  get z(): number {return this.coordinates[2] ? this.coordinates[2] : 0;}
  set x(val: number) {this.coordinates[0] = val;}
  set y(val: number) {this.coordinates[1] = val;}
  set z(val: number) {this.coordinates[2] = val;}
  get lng(): number {return this.coordinates[0];}
  get lat(): number {return this.coordinates[1];}
  get alt(): number {return this.coordinates[2] ? this.coordinates[2] : 0;}
  set lng(val: number) {this.coordinates[0] = val;}
  set lat(val: number) {this.coordinates[1] = val;}
  set alt(val: number) {this.coordinates[2] = val;}

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

  get cellType() {return SEMTYPEGIS.GISPOINT;}
  get defaultWidth() {return 200;}
  get defaultHeight() {return 100;}

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
    g.fillText('['+obj.x+';'+obj.y+']', x + 10, y + 10);
  }
}

//Point semantic type handler
export class GisPointHandler extends DG.ObjectHandler {
  get type() {return SEMTYPEGIS.GISPOINT;}

  isApplicable(obj: any) {return (obj instanceof GisPoint);}

  // getCanvasRenderer() {return null;}
  getCanvasRenderer() {return new GisPointCanvasRenderer();}
  getGridCellRenderer() {return new GisPointGridCellRenderer();}

  renderIcon() {return ui.iconFA('bullseye');}
  renderTooltip(obj: GisPoint) {return ui.divText(`[${obj.x} ; ${obj.y}]`);}
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
