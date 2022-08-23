/*
 * GIS VIEWER
 * TODO: add description here
*/
//base import
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//QJuery import
import $ from 'cash-dom';

//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer';
import {Coordinate} from '../src/gis-openlayer';
import {OLCallbackParam} from '../src/gis-openlayer';
import VectorLayer from 'ol/layer/Vector';


//type Pair<T, K> = [T, K]; //used Coordinate instead if it

export class GisViewer extends DG.JsViewer {
  latitudeColumnName: string;
  longitudeColumnName: string;
  valuesColumnName: string;
  renderType: string;
  weightedMarkers: boolean;
  markerSize: number;
  markerOpacity: number;
  heatmapRadius: number;
  heatmapBlur: number;

  isRowFocusing: boolean = true;

  //ui elements
  panelLeft: HTMLElement | null = null;
  panelRight: HTMLElement | null = null;
  panelTop: HTMLElement | null = null;
  panelBottom: HTMLElement | null = null;
  viewerContainer: HTMLElement | null = null;
  lblStatusCoords: HTMLElement | null = null;
  divLayersList: HTMLElement | null = null;
  leftPanelWidth: string = '110px';
  rightPanelWidth: string = '50px';

  initialized: boolean;
  ol: OpenLayers;
  layers = [];
  coordinates: Coordinate[] = [];
  values: Array<string | number> = [];

  constructor() {
    super();

    this.initialized = false;
    this.ol = new OpenLayers();

    // properties
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    this.valuesColumnName = this.string('valuesColumnName');
    this.weightedMarkers = this.bool('weightedMarkers', false);
    this.markerSize = this.int('markerSize', 5);
    this.markerOpacity = this.float('markerOpacity', 0.8, {editor: 'slider', min: 0.1, max: 1});

    this.heatmapRadius = this.int('heatmapRadius', 10);
    this.heatmapBlur = this.int('heatmapBlur', 20);
    // DG.Property
    //JsViewer

    this.renderType = this.string('renderType', 'markers');
    const renderTypeProp = this.getProperty('renderType');
    if (renderTypeProp) renderTypeProp.choices = ['markers', 'heat map'];
  }

  layerUIElement(num: number, layerName: string, layerId: string): HTMLElement {
    const l1 = ui.div(num+': '+layerName);
    l1.style.overflow = 'hidden';
    (l1.lastElementChild as HTMLElement).style.overflow = 'hidden';
    //TODO: add tooltip to base element with layer name

    function setupBtn(el: HTMLElement, layerName: string, layerId: string) {
    // el.style.border = '1px solid black';
      el.style.height = '14px';
      el.style.width = '18px';
      el.style.marginLeft = '4px';
      el.setAttribute('layerName', layerName);
      el.setAttribute('layerId', layerId);
    }
    //l1.style.border = '1px solid blue';
    // const btnVisible = ui.button(ui.iconFA('eye'), ()=>{
    const btnVisible = ui.iconFA('eye', (evt)=>{ //cogs
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (!layerId) return;
      const layer = this.ol.getLayerById(layerId);
      if (!layer) return;
      const isVisible = layer.getVisible();
      layer.setVisible(!isVisible);
      // if (!isVisible) divLayer.style.background = 'lightblue';
      // else divLayer.style.background = 'lightgray';
    }, 'Show/hide layer');
    setupBtn(btnVisible, layerName, layerId);
    //Setup button>>
    const btnSetup = ui.iconFA('cogs', (evt)=>{ //cogs
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (!layerId) return;
      const layer = this.ol.getLayerById(layerId) as VectorLayer<any>;
      if (!layer) return;
      const arrFeatures = this.ol.getFeaturesFromLayer(layer);
      if (arrFeatures) {
        const df = DG.DataFrame.fromObjects(arrFeatures);
        if (df) {
          df.name = layer.get('layerName');
          grok.shell.addTableView(df as DG.DataFrame);
        }
      }
    }, 'Layer properties');
    setupBtn(btnSetup, layerName, layerId);
    //Delete button>>
    const btnDelete = ui.iconFA('trash', (evt)=>{
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (layerId) this.ol.removeLayerById(layerId);
      this.updateLayersList();
    }, 'Delete layer');
    setupBtn(btnDelete, layerName, layerId);

    const panelButtons = ui.divH([btnVisible, btnSetup, btnDelete]);

    const divLayerUI = ui.divV([l1, panelButtons]);
    // divLayerUI.style.border = '1px solid darkgray';
    divLayerUI.style.overflow = 'hidden';
    return divLayerUI;
  }

  initUi(shortUI: boolean = false): HTMLElement {
    //create UI>>
    this.lblStatusCoords = ui.div('long, lat');
    this.lblStatusCoords.id = 'lbl-coord';
    const divStatusBody = ui.divH([ui.div('coord :'), this.lblStatusCoords]);
    this.panelBottom = ui.box(divStatusBody);
    this.panelBottom.style.maxHeight = '20px';
    // this.panelBottom.style.background = '#f2f2f5';
    this.panelBottom.style.border = 'solid 1px lightgray';

    this.panelTop = ui.box();
    this.panelTop.style.maxHeight = '36px';

    this.panelLeft = ui.box();
    this.panelLeft.style.maxWidth = '110px';
    this.panelLeft.style.minWidth = '60px';
    // this.panelLeft.style.border = 'solid 1px darkgray';

    //menu bar icons>>
    // const iconLayers = ui.iconImage('Layers', '/icons/layers-svgrepo.svg');
    //const leftPanelBtn = ui.button(iconLayers, ()=>{ //cogs
    const leftPanelBtn = ui.button(ui.iconFA('layer-group'), ()=>{ //cogs
      if (this.panelLeft) {
        if (this.panelLeft.style.visibility == 'visible') {
          this.panelLeft.style.visibility = 'hidden';
          this.leftPanelWidth = this.panelLeft.style.maxWidth;
          this.panelLeft.style.maxWidth = '2px';
          this.panelLeft.style.minWidth = '2px';
        } else {
          this.panelLeft.style.visibility = 'visible';
          this.panelLeft.style.maxWidth = this.leftPanelWidth;
          this.panelLeft.style.minWidth = '60px';
          this.updateLayersList();
        }
      }
      this.rootOnSizeChanged(this.root);
    });

    const btnRowFocusing = ui.button(ui.iconFA('bullseye'), ()=>{
      this.isRowFocusing = !this.isRowFocusing;
    });
    const heatmapBtn = ui.button('HM', ()=>{
      this.renderHeat();
      this.updateLayersList();
    });

    //add buttons to top menu panel
    this.panelTop.append(ui.divH([leftPanelBtn, heatmapBtn, btnRowFocusing]));
    //this.panelTop.append(ui.divH([leftPanelBtn, heatmapBtn]));

    //left panel icons>>
    this.divLayersList = ui.divV([]);
    this.panelLeft.append(this.divLayersList);

    leftPanelBtn.click();
    // this.panelLeft.style.visibility = 'hidden';
    // this.panelLeft.style.maxWidth = '2px';
    // this.panelLeft.style.minWidth = '2px';

    const body = ui.box();
    body.id = 'map-container';
    // body.style.border = 'solid 1px darkgray';
    //body.style.minWidth = '100px';

    this.viewerContainer = ui.splitV(
      [this.panelTop,
        ui.splitH([this.panelLeft, body], null, true),
        this.panelBottom]);

    // this.viewerContainer.style.border = 'solid 2px lightgray';
    this.root.appendChild(this.viewerContainer);

    return this.viewerContainer;
  }

  init() {
    try {
      ui.setUpdateIndicator(this.root, true);
      this.initUi();
      // ui.setUpdateIndicator(this.panelTop!, true);
      this.ol.initMap('map-container');

      //TODO: refactor here>
      this.ol.markerSize = this.markerSize;
      this.ol.markerOpacity = this.markerOpacity;
      this.ol.weightedMarkers = this.weightedMarkers;
      this.ol.heatmapRadius = this.heatmapRadius;
      this.ol.heatmapBlur = this.heatmapBlur;

      this.updateLayersList();

      //subscribe to events
      this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
      this.subs.push(ui.onSizeChanged((this.panelLeft as HTMLElement)).subscribe(this.rootOnSizeChanged.bind(this)));
      //setup callbacks
      this.ol.setMapPointermoveCallback(this.showCoordsInStatus.bind(this));
      // this.ol.setMapClickCallback(this.showCoordsInStatus.bind(this));
      this.ol.setMapClickCallback(this.handlerOnMapClick.bind(this));

      this.initialized = true;
    }
    catch (e: any) {
      this.initialized = false;
      grok.shell.error(e.toString());
      this.root.appendChild(
        ui.divV([ui.div('Error loading GIS map!'), ui.div(e.toString())]));
    }
    finally {
      ui.setUpdateIndicator(this.root, false);
      // ui.setUpdateIndicator(this.panelTop!, true);
    }
  }

  updateLayersList() {
    if ((!this.ol) || (!this.panelLeft) || (!this.divLayersList)) return;
    // let htmlStyle: DG.ElementOptions = { };

    while (this.divLayersList.lastChild)
      this.divLayersList.removeChild(this.divLayersList.lastChild);

    const layersArr = this.ol.getLayersList();

    for (let i = 0; i < layersArr.length; i++) {
      let layerName = layersArr[i].get('layerName');
      if ((layerName === undefined) || (layerName === null) ) layerName = '';
      const layerId = layersArr[i].get('layerId');

      // const divLayer = ui.div(i+' '+layerName, {style: {'border': 'solid 1px lightgray'}});
      const divLayer = this.layerUIElement(i, layerName, layerId);
      const isVisible = layersArr[i].getVisible();
      divLayer.setAttribute('layerName', layerName);
      divLayer.setAttribute('layerId', layerId);
      // if (isVisible) divLayer.style.background = 'lightblue';
      // else divLayer.style.background = 'lightgray';

      divLayer.onclick = (evt)=>{
        // const divLayer = (evt.currentTarget as HTMLElement);
        // // const layerName = divLayer.getAttribute('layerName');
        // // if (!layerName) return;
        // // const layer = this.ol.getLayerByName(layerName);
        // const layerId = divLayer.getAttribute('layerId');
        // if (!layerId) return;
        // const layer = this.ol.getLayerById(layerId);
        // if (!layer) return;
        // const isVisible = layer.getVisible();
        // layer.setVisible(!isVisible);
        // if (!isVisible) divLayer.style.background = 'lightblue';
        // else divLayer.style.background = 'lightgray';
      };

      this.divLayersList.append(divLayer);
      // this.divLayersList.append(this.layerUIElement(layerName, layerId));
    }
  }

  showCoordsInStatus(p: OLCallbackParam): void {
    if (this.lblStatusCoords) {
      if (p)
        this.lblStatusCoords.innerText = (p.coord[0]).toFixed(3) + ' / ' + (p.coord[1]).toFixed(3);
    }
  }

  handlerOnMapClick(p: OLCallbackParam): void {
    if (this.lblStatusCoords) {
      if (p) {
        if (this.isRowFocusing) {
          // this.lblStatusCoords.innerText = (p.coord[0]).toFixed(3) + ' / ' + (p.coord[1]).toFixed(3);
          this.dataFrame.rows.select(
            (row) => ((row[this.latitudeColumnName] == p.coord[0]) && (row[this.longitudeColumnName] == p.coord[1])));
        }
      }
    }
  }

  private rootOnSizeChanged(args: any): void {
    this.ol.olMap.setSize([this.root.clientWidth, this.root.clientHeight]);
    //this.ol.olMap.updateSize();
    setTimeout( function(m) {m.updateSize();}, 200, this.ol.olMap);
  }

  onTableAttached(): void {
    this.init();

    const latLngColumns = this.dataFrame.columns.bySemTypesExact([DG.SEMTYPE.LATITUDE, DG.SEMTYPE.LONGITUDE]);
    if (latLngColumns != null && this.latitudeColumnName == null && this.longitudeColumnName == null) {
      this.latitudeColumnName = latLngColumns[0].name;
      this.longitudeColumnName = latLngColumns[1].name;
    }

    // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    // this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.render(true);
  }

  onPropertyChanged(prop: DG.Property): void {
    if (!this.initialized) return;
    //property.name

    this.render();
  }

  detach(): void {
    //TODO - release map
    // this.map.remove();
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  render(fit: boolean = false): void {
    //TODO: refactor this>
    this.ol.markerSize = this.markerSize;
    this.ol.markerOpacity = this.markerOpacity;
    this.ol.weightedMarkers = this.weightedMarkers;
    this.ol.heatmapRadius = this.heatmapRadius;
    this.ol.heatmapBlur = this.heatmapBlur;

    if (this.latitudeColumnName == null || this.longitudeColumnName == null)
      return;

    this.getCoordinates(); //not only coordinates but a data at all

    if (this.renderType === 'heat map')
      this.renderHeat();
    else if (this.renderType === 'markers')
      this.renderMarkers();

    // if (fit)
    //   this.map.fitBounds(this.coordinates);
    this.updateLayersList();
  }

  getCoordinates(): void {
    this.coordinates.length = 0;
    this.values.length = 0;
    const indexes = this.dataFrame.filter.getSelectedIndexes();
    let val: Int32Array | Float32Array | Float64Array | Uint32Array;

    const lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    const lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();
    // val = this.dataFrame.getCol(this.valuesColumnName).getRawData();

    if ((!lat) || (!lon)) return;

    const col = this.dataFrame.getCol(this.valuesColumnName);
    //TODO: change it to filling array of objects with all table data
    if (col) val = col.getRawData();
    else {
      val = new Float32Array(lat.length);
      val.fill(1);
    }

    for (let i = 0; i < indexes.length; i++) {
      if ((i < lat.length) && (i < lon.length)) {
        this.coordinates.push([lon[indexes[i]], lat[indexes[i]]]);
        this.values.push(val[indexes[i]]);
      }
    }
  }

  renderHeat(): void {
    this.getCoordinates();
    let colName = this.valuesColumnName;
    if (colName === 'null') colName = this.dataFrame.name;
    if (colName === null) colName = this.dataFrame.name;

    let layer = this.ol.getLayerByName('HL: ' + colName);
    if (!layer) layer = this.ol.addNewHeatMap('HL: ' + colName);

    this.ol.clearLayer(layer);
    for (let i = 0; i < this.coordinates.length; i++)
      this.ol.addPoint(this.coordinates[i], this.values[i], layer);
  }

  renderMarkers(): void {
    this.ol.clearLayer(); //TODO: clear exact layer
    for (let i = 0; i < this.coordinates.length; i++)
      this.ol.addPoint(this.coordinates[i], this.values[i]);
  }
}
