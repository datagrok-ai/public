/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const SEMTYPEGIS = {
  LONGITUDE: 'gis-longitude',
  LATIITUDE: 'gis-latitude',
  ALTITUDE: 'gis-altitude',
  GISPOINT: 'gis-point',
  GISAREA: 'gis-gis_area',
};

export type gisCoordinate = [number, number, number?];
export type gisFeatureProperties = {[x: string]: any};

export class GisArea {
  semtype: string = SEMTYPEGIS.GISAREA;
  coordinates: Array<gisCoordinate>; //Array<[number, number, number?]>
  attributes: gisFeatureProperties = {};
  //TODO: add x, y coordinate-getters that incapsulate lng, lat
  constructor(coord: Array<gisCoordinate>, attr?: gisFeatureProperties) {
    this.coordinates = coord;
    if (attr) this.attributes = attr;
  }
}

export class GisAreaCanvasRenderer extends DG.CanvasRenderer {
  get defaultWidth(): number | null {
    return null;
  }
  get defaultHeight(): number | null {
    return null;
  }
  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    obj: GisArea, context: any): void {
    if (obj.coordinates.length == 0) return;
    g.fillStyle = '#FF0000';
    for (let i = 0; i < obj.coordinates.length-1; i++) {
      g.moveTo(obj.coordinates[i][0], obj.coordinates[i][1]);
      g.lineTo(obj.coordinates[i+1][0], obj.coordinates[i+1][1]);
      g.stroke();
      // panelProperties.appendChild(ui.divText(`[${obj.coordinates[i][0]} ; ${obj.coordinates[i][1]}]`));
      // g.fillText('['+obj.coordinates[i][0]+';'+obj.coordinates[i][1]+']', x + 10, y + 10);
    }
    //g.fillText('['+obj.longitude+';'+obj.latitude+']', x + 10, y + 10);
  }
}

//Point semantic type handler
export class GisAreaHandler extends DG.ObjectHandler {
  get type() {return SEMTYPEGIS.GISPOINT;}

  isApplicable(obj: any) {return (obj instanceof GisArea);}

  // getCanvasRenderer() {return null;}
  getCanvasRenderer() {return new GisAreaCanvasRenderer();}
  // getGridCellRenderer() {return new GisAreaGridCellRenderer();}

  renderIcon() {return ui.iconFA('bullseye');}
  renderTooltip(obj: GisArea) {return ui.divText(`Area of [${obj.coordinates.length}] vertexes`);}
  renderProperties(obj: GisArea) {
    const panelProperties = ui.divV([]);
    panelProperties.appendChild(ui.divText(`Properties for Area`));
    for (let i = 0; i < obj.coordinates.length; i++)
      panelProperties.appendChild(ui.divText(`[${obj.coordinates[i][0]} ; ${obj.coordinates[i][1]}]`));
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
  semtype: string = 'gis_point';
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
}

export class GisPointGridCellRenderer extends DG.GridCellRenderer {
  get cellType() {return SEMTYPEGIS.GISPOINT;}

  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    g.fillStyle = 'darkgray';
    // gridCell.ob
    // g.fillText('['+this.longitude+';'+this.latitude+']', x + 5, y + 5);
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
