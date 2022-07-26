/*
  GIS VIEWER
  TODO: add description here
*/
//base import
//import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer';

type Pair<T, K> = [T, K];

export class GisViewer extends DG.JsViewer {
  latitudeColumnName: string;
  longitudeColumnName: string;
  // renderType = this.string('renderType', 'heat map');
  // markerSize = this.int('markerSize', 10);
  // markerOpacity = this.float('markerOpacity', 0.8);
  //getProperty('renderType').choices = ['markers', 'heat map'];

  initialized: boolean;
  ol: OpenLayers;
  layers = [];
  coordinates: Pair<number, number>[] = [];

  constructor() {
    super();

    this.initialized = false;
    this.ol = new OpenLayers();

    // properties
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    // this.renderType = this.string('renderType', 'heat map');
    // this.markerSize = this.int('markerSize', 10);
    // this.markerOpacity = this.float('markerOpacity', 0.8);
    // this.getProperty('renderType').choices = ['markers', 'heat map'];

    // this.layers = [];
    // this.coordinates = [];

    // ui.onSizeChanged(this.root).subscribe((_) => this.map.invalidateSize());
  }

  initUi() {
    let htmlStyle: DG.ElementOptions = { };
    htmlStyle = {style: {'border': 'solid 1px darkgray'}};
    const panelTop = ui.divH([ui.div('MENU')], htmlStyle);
    const panelLeft = ui.divV([], htmlStyle);
    const panelRight = ui.divV([], htmlStyle);
    const panelBottom = ui.divH([ui.div('status bar')], htmlStyle);
    // panelTop.setAttribute('id', 'gis-panel-top'); //TODO: maybe we need to use id further
    panelBottom.style.background = '#f2f2f5';

    //const boxMap = ui.box(null, 'd4-viewer-host');
    htmlStyle = {style: {'width': '100%', 'height': '100%', 'border': 'solid 1px red'}};
    //const boxMap = ui.box(null, htmlStyle);
    //const boxMap = ui.div([], 'd4-viewer-host');
    const boxMap = ui.div([], htmlStyle);
    boxMap.id = 'map-container';
    //htmlStyle = {style: {'border': 'solid 1px blue'}};
    htmlStyle = {style: {'width': '100%', 'height': '100%', 'border': 'solid 1px blue'}};
    const boxCenter = ui.divH([panelLeft, boxMap, panelRight], htmlStyle);
    //const boxCenter = ui.divH([panelLeft, boxMap, panelRight], 'd4-viewer-host');

    const divViewerContainer = ui.divV([panelTop, boxCenter, panelBottom], 'd4-viewer-host');

    this.root.appendChild(divViewerContainer);
  }

  init() {
    this.initUi();
    this.ol.initMap('map-container');

    this.initialized = true;
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
    // this.map.remove();
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  render(fit: boolean = false) {
    // for (const layer of this.layers)
    //   this.map.removeLayer(layer);
    // this.layers.length = 0;

    if (this.latitudeColumnName == null || this.longitudeColumnName == null)
      return;

    this.getCoordinates();
    this.renderMarkers();

    // if (this.renderType === 'heat map')
    //   this.renderHeat();
    // else if (this.renderType === 'markers')
    //   this.renderMarkers();

    // if (fit)
    //   this.map.fitBounds(this.coordinates);
  }

  getCoordinates() {
    this.coordinates.length = 0;
    const indexes = this.dataFrame.filter.getSelectedIndexes();
    const lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    const lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();

    for (let i = 0; i < indexes.length; i++) {
      if ((i < lat.length ) && (i < lon.length ))
        this.coordinates.push([lon[indexes[i]], lat[indexes[i]]]);
    }
  }

  renderHeat() {
    // this.layers.push(L.heatLayer(this.coordinates, {radius: this.markerSize}).addTo(this.map));
  }

  renderMarkers() {
    for (let i = 0; i < this.coordinates.length; i++)
      this.ol.addPoint(this.coordinates[i]);

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
