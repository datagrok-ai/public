/*
 * GIS VIEWER
 * TODO: add description here
*/
//base import
//import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//QJuery import
import $ from 'cash-dom';

//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer';
import {Coordinate} from '../src/gis-openlayer';
import {OLCallbackParam} from '../src/gis-openlayer';


type Pair<T, K> = [T, K];

export class GisViewer extends DG.JsViewer {
  latitudeColumnName: string;
  longitudeColumnName: string;
  valuesColumnName: string;
  renderType: string;
  markerSize: number;
  markerOpacity: number;
  weightedMarkers: boolean;

  //ui elements
  panelLeft: HTMLElement | null = null;
  panelRight: HTMLElement | null = null;
  panelTop: HTMLElement | null = null;
  panelBottom: HTMLElement | null = null;
  viewerContainer: HTMLElement | null = null;
  lblStatusCoords: HTMLElement | null = null;
  divLayersList: HTMLElement | null = null;

  initialized: boolean;
  ol: OpenLayers;
  layers = [];
  coordinates: Coordinate[] = [];
  values: Array<string | number> = [];
  // coordinates: Pair<number, number>[] = [];

  constructor() {
    super();

    this.initialized = false;
    this.ol = new OpenLayers();

    // properties
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    this.valuesColumnName = this.string('valuesColumnName');
    this.markerSize = this.int('markerSize', 10);
    this.markerOpacity = this.float('markerOpacity', 0.8);
    this.weightedMarkers = this.bool('weightedMarkers', true);

    this.renderType = this.string('renderType', 'markers');
    const renderTypeProp = this.getProperty('renderType');
    if (renderTypeProp) renderTypeProp.choices = ['markers', 'heat map'];
  }

  initUi(shortUI: boolean = false): HTMLElement {
    //create UI>>
    this.lblStatusCoords = ui.div('long, lat');
    this.lblStatusCoords.id = 'lbl-coord';
    const divStatusBody = ui.divH([ui.div('status bar: '), this.lblStatusCoords]);
    this.panelBottom = ui.box(divStatusBody);
    this.panelBottom.style.maxHeight = '20px';
    this.panelBottom.style.background = '#f2f2f5';
    this.panelBottom.style.border = 'solid 1px darkgray';

    this.panelTop = ui.box();
    this.panelTop.style.maxHeight = '36px';

    this.panelLeft = ui.box();
    this.panelLeft.style.maxWidth = '90px';
    this.panelLeft.style.border = 'solid 1px darkgray';

    this.panelRight = ui.box();
    this.panelRight.style.maxWidth = '10px';
    //this.panelRight.style.border = 'solid 1px darkgray';
    //menu bar icons>>
    const leftPanelBtn = ui.button('Left panel', ()=>{
      // this.panelLeft.style.visibility = (this.panelLeft.style.visibility == 'visible') ? 'hidden' : 'visible';
      if (this.panelLeft) {
        if (this.panelLeft.style.visibility == 'visible') {
          this.panelLeft.style.visibility = 'hidden';
          this.panelLeft.style.maxWidth = '2px';
        } else {
          this.panelLeft.style.visibility = 'visible';
          this.panelLeft.style.maxWidth = '90px';
        }
      }
      this.rootOnSizeChanged(null);
      //$(this.panelLeft).hide();
    });
    const rightPanelBtn = ui.button('Right panel', ()=>{
      // if (this.panelRight)
      //   this.panelRight.style.visibility = (this.panelRight.style.visibility == 'visible') ? 'hidden' : 'visible';
      if (this.panelRight) {
        if (this.panelRight.style.visibility == 'visible') {
          this.panelRight.style.visibility = 'hidden';
          this.panelRight.style.maxWidth = '2px';
        } else {
          this.panelRight.style.visibility = 'visible';
          this.panelRight.style.maxWidth = '10px';
        }
      }
      this.rootOnSizeChanged(null);
      //$(this.panelRight).hide();
    });
    const heatmapBtn = ui.button('HM', ()=>{
      this.renderHeat();
      //$(this.panelRight).hide();
    });
    //add buttont to top menu panel
    this.panelTop.append(ui.divH([leftPanelBtn, rightPanelBtn, heatmapBtn]));

    //left panel icons>>
    const l1 = ui.div(' <l1>');
    const l2 = ui.div(' <l2>');
    this.divLayersList = ui.divV([l1, l2]);
    this.panelLeft.append(this.divLayersList);

    this.panelLeft.style.visibility = 'visible';
    this.panelRight.style.visibility = 'visible';
    // $(this.panelLeft).hide();
    // $(this.panelRight).hide();

    const body = ui.box();
    body.id = 'map-container';

    this.viewerContainer = ui.splitV(
      [this.panelTop,
        ui.splitH([this.panelLeft, body, this.panelRight]),
        this.panelBottom]);

    this.viewerContainer.style.border = 'solid 2px lightgray';
    this.root.appendChild(this.viewerContainer);

    return this.viewerContainer;
  }

  init() {
    this.initUi();
    this.ol.initMap('map-container');

    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

    this.ol.setMapPointermoveCallback(this.showCoordsInStatus);
    this.ol.labelStatus = this.lblStatusCoords; //TODO: remove this construction - all interactions through callback!!

    this.updateLayersList();

    this.initialized = true;
  }

  updateLayersList() {
    if (!this.ol) return;
    if (!this.panelLeft) return;
    if (!this.divLayersList) return;

    //while (this.divLayersList.lastChild) this.divLayersList.lastChild.remove;
    while (this.divLayersList.lastChild) this.divLayersList.removeChild(this.divLayersList.lastChild);
    const layersNames = this.ol.getLayersList();

    for (let i = 0; i < layersNames.length; i++) {
      // const l2 = ui.div(layersNames[i]);
      this.divLayersList.append(ui.div(layersNames[i], {style: {'border': 'solid 1px blue'}}));
    }
  }

  showCoordsInStatus(p: OLCallbackParam): void {
    if (this.lblStatusCoords) {
      if (p)
        this.lblStatusCoords.innerText = p.coord[0] + ', ' + p.coord[1];
    }
  }

  private rootOnSizeChanged(args: any): void {
    //this.ol.olMap.setSize(this.root.clientWidth, this.root.clientHeight);
    //this.ol.olMap.updateSize();
    setTimeout( function(m) {m.updateSize();}, 200, this.ol.olMap);
  }

  onTableAttached() {
    this.init();

    const latLngColumns = this.dataFrame.columns.bySemTypesExact([DG.SEMTYPE.LATITUDE, DG.SEMTYPE.LONGITUDE]);
    if (latLngColumns != null && this.latitudeColumnName == null && this.longitudeColumnName == null) {
      this.latitudeColumnName = latLngColumns[0].name;
      this.longitudeColumnName = latLngColumns[1].name;
    }

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.render(true);
  }

  onPropertyChanged(prop: any) {
    if (this.initialized)
      this.render();
  }

  detach() {
    //TODO - release map
    // this.map.remove();
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  render(fit: boolean = false) {
    // for (const layer of this.layers)
    //   this.map.removeLayer(layer);
    // this.layers.length = 0;

    if (this.latitudeColumnName == null || this.longitudeColumnName == null)
      return;

    this.getCoordinates(); //not only coordinates but a data at all
    // this.renderMarkers();

    if (this.renderType === 'heat map')
      this.renderHeat();
    else if (this.renderType === 'markers')
      this.renderMarkers();

    // if (fit)
    //   this.map.fitBounds(this.coordinates);
    this.updateLayersList();
  }

  getCoordinates() {
    this.coordinates.length = 0;
    this.values.length = 0;
    const indexes = this.dataFrame.filter.getSelectedIndexes();
    let lat: Int32Array | Float32Array | Float64Array | Uint32Array;
    let lon: Int32Array | Float32Array | Float64Array | Uint32Array;
    let val: Int32Array | Float32Array | Float64Array | Uint32Array;

    lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();
    // val = this.dataFrame.getCol(this.valuesColumnName).getRawData();

    // col = this.dataFrame.getCol(this.latitudeColumnName);
    // if (col) lat = col.getRawData();
    // col = this.dataFrame.getCol(this.longitudeColumnName);
    // if (col) lon = col.getRawData();
    const col = this.dataFrame.getCol(this.valuesColumnName);
    //TODO: change it to filling array of objects with all table data
    if (col) val = col.getRawData();
    else {
      val = new Int32Array(lat.length); //lat.copyWithin(0, 0);
      val.fill(1);
    }

    for (let i = 0; i < indexes.length; i++) {
      if ((i < lat.length ) && (i < lon.length )) {
        this.coordinates.push([lon[indexes[i]], lat[indexes[i]]]);
        this.values.push(val[indexes[i]]);
      }
    }
  }

  renderHeat() {
    // this.layers.push(L.heatLayer(this.coordinates, {radius: this.markerSize}).addTo(this.map));
    this.getCoordinates();
    const layer = this.ol.addNewHeatMap('heatLayer');
    for (let i = 0; i < this.coordinates.length; i++)
      this.ol.addPoint(this.coordinates[i], this.values[i], layer);
  }

  renderMarkers() {
    this.ol.clearLayer(); //TODO: clear exact layer
    for (let i = 0; i < this.coordinates.length; i++)
      this.ol.addPoint(this.coordinates[i], this.values[i]);

    // let markerOptions = {
    //   radius: this.markerSize,
    //   fillColor: "#ff7800",
    //   color: "#000",
    //   weight: 1,
    //   opacity: 1,
    //   fillOpacity: this.markerOpacity
    // };

    // let markers = [];
    // for (let i = 0; i < this.coordinates.length; i++)
    //   markers.push(L.circleMarker(this.coordinates[i], markerOptions));
    // this.layers.push(L.featureGroup(markers).addTo(this.map));
  }
}
