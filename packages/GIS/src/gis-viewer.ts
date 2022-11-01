/*
 * GIS VIEWER
 * this is the viewer module which manage viewer UI, map based on open-layers class
*/
//base import
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//QJuery import
import $ from 'cash-dom';

//GIS semantic types import
import {SEMTYPEGIS} from '../src/gis-semtypes';

//OpenLayers functionality import
import {OpenLayers} from '../src/gis-openlayer';
import {Coordinate} from '../src/gis-openlayer';
import {OLCallbackParam} from '../src/gis-openlayer';
import VectorLayer from 'ol/layer/Vector';
import * as OLProj from 'ol/proj';
import Feature from 'ol/Feature';
import { Point } from 'ol/geom';
import { COLOR_CODING_TYPE } from 'datagrok-api/dg';
import { forEach } from 'jszip';
// import { numberSafeCompareFunction } from 'ol/array';

function toStringColor(num : number, opacity?: number) : string {
  num >>>= 0;
  const b = num & 0xFF;
  const g = (num & 0xFF00) >>> 8;
  const r = (num & 0xFF0000) >>> 16;
  const a = opacity ? opacity : 1;
  return 'rgba(' + [r, g, b, a].join(',') + ')';
}

export class GisViewer extends DG.JsViewer {
  currentLayer: string;
  latitudeColumnName: string;
  longitudeColumnName: string;
  sizeColumnName: string;
  colorColumnName: string;
  labelsColumnName: string;
  markersColumnName: string;
  markerOpacity: number;
  defaultColor: number;
  selectedColor: number;
  markerMinColor: number;
  markerMaxColor: number;
  markerDefaultSize: number;
  markerMinSize: number;
  markerMaxSize: number;
  // gradientColoring: boolean;
  // gradientSizing: boolean;
  renderType: string;
  heatmapRadius: number;
  heatmapBlur: number;
  showTooltip: boolean;


  isRowFocusing: boolean = true;
  isShortUI: boolean = true;

  //ui elements
  panelLeft: HTMLElement | null = null;
  leftPanelBtn: HTMLElement | null = null;
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
  features: Array<Feature> = [];
  coordinates: Coordinate[] = [];
  labelValues: Array<string | number> = [];
  colorValues: Array<number> = [];
  sizeValues: Array<number> = [];
  indexValues: Array<number> = [];

  constructor() {
    super();

    this.initialized = false;
    this.ol = new OpenLayers(this);

    // properties
    this.currentLayer = this.string('currentLayer', 'BaseLayer', {category: 'Layers', choices: ['l1','l2','l3','l4']});
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    this.colorColumnName = this.string('colorColumnName');
    this.sizeColumnName = this.string('sizeColumnName');
    this.labelsColumnName = this.string('labelsColumnName', '', {userEditable: false});
    this.markersColumnName = this.string('markersColumnName', '', {userEditable: false});

    this.markerOpacity = this.float('markerOpacity', 80, {category: 'Markers', editor: 'slider', min: 0, max: 100});
    // this.gradientSizing = this.bool('gradientSizing', false, {category: 'Markers'});
    this.markerDefaultSize = this.int('markerDefaultSize', 5, {category: 'Markers', min: 1, max: 30});
    //DONE: there is no need in valitators: Min can be > Max
    this.markerMinSize = this.int('markerMinSize', 2, {category: 'Markers', min: 1, max: 30});
    this.markerMaxSize = this.int('markerMaxSize', 15, {category: 'Markers', min: 1, max: 30});

    // this.gradientColoring = this.bool('gradientColoring', false, {category: 'Markers'});
    this.defaultColor = this.int('defaultColor', 0x1f77b4, {category: 'Markers'});
    this.selectedColor = this.int('selectedColor', 0xff8c00, {category: 'Markers'});
    this.markerMinColor = this.int('markerMinColor', 0x0000ff, {category: 'Markers', userEditable: false});
    this.markerMaxColor = this.int('markerMaxColor', 0xff0000, {category: 'Markers', userEditable: false}); //hexToRGB

    this.heatmapRadius = this.int('heatmapRadius', 5,
      {category: 'Heatmap', description: 'Heatmap radius', min: 1, max: 10, userEditable: true});
    this.heatmapBlur = this.int('heatmapBlur', 20,
      {category: 'Heatmap', description: 'Heatmap radius', min: 1, max: 50, userEditable: true});
    // DG.Property
    //JsViewer


    this.showTooltip = this.bool('showTooltip', false);

    this.renderType = this.string('renderType', 'markers');
    const renderTypeProp = this.getProperty('renderType');
    if (renderTypeProp)
      renderTypeProp.choices = ['markers', 'heat map'];
  }

  layerUIElement(num: number, layerName: string, layerId: string, isVisible: boolean): HTMLElement {
    const l1 = ui.div(num+': '+layerName);
    l1.style.overflow = 'hidden';
    l1.style.flexWrap = 'nowrap';
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
    const btnVisible = ui.iconFA(isVisible ? 'eye' : 'eye-slash', (evt)=>{ //cogs
      evt.stopImmediatePropagation();
      // evt.stopPropagation();
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (!layerId)
        return;
      const layer = this.ol.getLayerById(layerId);
      if (!layer)
        return;
      const isVisible = layer.getVisible();
      layer.setVisible(!isVisible);
      this.updateLayersList();
      // if (!isVisible) divLayer.style.background = 'lightblue';
      // else divLayer.style.background = 'lightgray';
    }, 'Show/hide layer');
    //temporary this is the save layer button
    setupBtn(btnVisible, layerName, layerId);
    //Setup button>>
    const btnSetup = ui.iconFA('download', (evt)=>{ //cogs
      evt.stopImmediatePropagation();
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (!layerId)
        return;
      const layer = this.ol.getLayerById(layerId) as VectorLayer<any>;
      if (!layer)
        return;
      const arrPreparedToDF: any[] = [];
      const arrFeatures = this.ol.getFeaturesFromLayer(layer);
      if (arrFeatures) {
        for (let i = 0; i < arrFeatures.length; i++) {
          // const gisObj = OpenLayers.gisObjFromGeometry(arrFeatures[i]);
          const newObj = arrFeatures[i].getProperties();
          if (newObj.hasOwnProperty('geometry'))
            delete newObj.geometry;
          if (arrFeatures[i].getId())
            newObj.id_ = arrFeatures[i].getId();
          newObj.gisObject = OpenLayers.gisObjFromGeometry(arrFeatures[i]);

          arrPreparedToDF.push(newObj);
        }
        // const df = DG.DataFrame.fromObjects(arrFeatures);
        const df = DG.DataFrame.fromObjects(arrPreparedToDF);
        if (df) {
          const gisCol = df.col('gisObject');
          if (gisCol)
            gisCol.semType = SEMTYPEGIS.GISAREA; //SEMTYPEGIS.GISOBJECT;
          df.name = layer.get('layerName');

          // df.toCsv()
          // this.view.addTableView(df as DG.DataFrame);
          grok.shell.addTableView(df as DG.DataFrame);

          // this.view.addTableView(df as DG.DataFrame);
        }
      }
    }, 'Layer properties');
    setupBtn(btnSetup, layerName, layerId);
    //Delete button>>
    const btnDelete = ui.iconFA('trash-alt', (evt)=>{ //trash
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (layerId)
        this.ol.removeLayerById(layerId);
      this.updateLayersList();
    }, 'Delete layer');
    setupBtn(btnDelete, layerName, layerId);

    const panelButtons = ui.divH([btnVisible, btnSetup, btnDelete]);

    const divLayerUI = ui.divV([l1, panelButtons]);
    if (layerName == this.currentLayer)
      divLayerUI.style.border = '2px solid blue';
    else
      divLayerUI.style.border = '1px solid lightgray';
    divLayerUI.style.overflow = 'hidden';
    return divLayerUI;
  }

  switchUI(shortUI: boolean) {
    if (this.panelTop) {
      this.panelTop.style.visibility = shortUI ? 'hidden' : 'visible';
      this.panelTop.style.maxHeight = shortUI ? '1px' : '30px';
    }
    if (this.panelBottom) {
      this.panelBottom.style.maxHeight = shortUI ? '2px' : '2px'; //temporary completely hide bottom panel
      this.panelBottom.style.visibility = shortUI ? 'hidden' : 'hidden'; //temporary completely hide bottom panel
      // this.panelBottom.style.visibility = shortUI ? 'hidden' : 'visible';
      // this.panelBottom.style.maxHeight = shortUI ? '2px' : '20px';
    }
    if ((this.panelLeft) && (this.leftPanelBtn)) {
      this.panelLeft.style.visibility = shortUI ? 'visible' : 'hidden';
      this.leftPanelBtn.click();
    }
  }

  initUi(shortUI: boolean = true): HTMLElement {
    this.isShortUI = shortUI;
    //create UI>>
    this.lblStatusCoords = ui.div('long, lat');
    this.lblStatusCoords.id = 'lbl-coord';
    const divStatusBody = ui.divH([ui.div('coord :'), this.lblStatusCoords]);
    this.panelBottom = ui.box(divStatusBody);
    this.panelBottom.style.maxHeight = shortUI ? '2px' : '2px'; //temporary hide bottom panel
    this.panelBottom.style.visibility = shortUI ? 'hidden' : 'hidden'; //temporary hide bottom panel

    this.panelTop = ui.box();
    this.panelTop.style.maxHeight = shortUI ? '1px' : '30px';
    this.panelTop.style.visibility = shortUI ? 'hidden' : 'visible';

    this.panelLeft = ui.box();
    this.panelLeft.style.maxWidth = '110px';
    this.panelLeft.style.minWidth = '70px';

    //menu bar icons>>
    this.leftPanelBtn = ui.button(ui.iconFA('layer-group'), ()=>{ //cogs
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
    this.leftPanelBtn.style.margin = '1px';

    const btnRowFocusing = ui.button(ui.iconFA('bullseye'), ()=>{
      this.isRowFocusing = !this.isRowFocusing;
      // this.dataFrame.
    });
    btnRowFocusing.style.margin = '1px';
    const heatmapBtn = ui.button('HM', ()=>{
      this.getCoordinates();
      this.renderHeat(this.features);
      this.updateLayersList();
    }, 'Build heat-map for data');
    heatmapBtn.style.margin = '1px';

    //add buttons to top menu panel
    // this.panelTop.append(ui.divH([leftPanelBtn, heatmapBtn, btnRowFocusing]));
    this.panelTop.append(ui.divH([this.leftPanelBtn, heatmapBtn]));

    //left panel icons>>
    this.divLayersList = ui.divV([]);
    this.panelLeft.append(this.divLayersList);

    this.panelLeft.style.visibility = 'visible';
    this.leftPanelBtn.click();

    const body = ui.box();
    body.id = 'map-container';
    body.style.maxWidth = '100%';
    body.style.maxHeight = '100%';

    this.viewerContainer = ui.splitV(
      [this.panelTop,
        ui.splitH([this.panelLeft, body], null, true),
        this.panelBottom]);

    this.root.appendChild(this.viewerContainer);

    //setup context menu
    this.onContextMenu.subscribe((menu) => {
      if (this.isShortUI == true) {
        menu.item('Extended UI', () => {
          this.isShortUI = false;
          this.switchUI(this.isShortUI);
        });
      } else {
        menu.item('Shortened UI', () => {
          this.isShortUI = true;
          this.switchUI(this.isShortUI);
        });
      }
    });

    return this.viewerContainer;
  }

  init() {
    try {
      ui.setUpdateIndicator(this.root, true);
      this.ol.useWebGL = true; //choose true if we want to use WebGL renderer

      this.initUi();
      // ui.setUpdateIndicator(this.panelTop!, true);
      this.ol.initMap('map-container');

      this.updateOpenLayerProperties();

      this.updateLayersList();

      //subscribe to events
      this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
      this.subs.push(ui.onSizeChanged((this.panelLeft as HTMLElement)).subscribe(this.rootOnSizeChanged.bind(this)));
      //setup callbacks
      this.ol.setMapPointermoveCallback(this.showMarkerTooltip.bind(this));
      // this.ol.setMapClickCallback(this.showCoordsInStatus.bind(this));

      this.ol.setMapSelectionCallback(this.handlerOnMarkerSelect.bind(this));
      this.ol.setMapClickCallback(this.handlerOnMapClick.bind(this));
      this.ol.setMapRefreshCallback(this.updateLayersList);

      this.initialized = true;
    } catch (e: any) {
      this.initialized = false;
      grok.shell.error(e.toString());
      this.root.appendChild(
        ui.divV([ui.div('Error loading GIS map!'), ui.div(e.toString())]));
    } finally {
      ui.setUpdateIndicator(this.root, false);
      setTimeout(() => {
        grok.shell.o = this;
        grok.shell.windows.showProperties = true;
      }, 500);
    }
  }

  updateLayersList(): void {
    if ((!this.ol) || (!this.panelLeft) || (!this.divLayersList))
      return;

    while (this.divLayersList.lastChild)
      this.divLayersList.removeChild(this.divLayersList.lastChild);

    const layersArr = this.ol.getLayersList();
    const arrLayerNames: [string] = [''];
    arrLayerNames.pop();

    for (let i = 0; i < layersArr.length; i++) {
      let layerName = layersArr[i].get('layerName');
      if ((layerName === undefined) || (layerName === null))
        layerName = '';
      const layerId = layersArr[i].get('layerId');
      const isVisible = layersArr[i].getVisible();

      arrLayerNames.push(layerName);
      const divLayer = this.layerUIElement(i, layerName, layerId, isVisible);
      divLayer.setAttribute('layerName', layerName);
      divLayer.setAttribute('layerId', layerId);

      divLayer.onclick = (evt)=>{
        const layerName = divLayer.getAttribute('layerName');
        const layerID = divLayer.getAttribute('layerId');
        if (!layerName || !layerID)
          return;
        this.currentLayer = layerName;
        this.setOptions({currentLayer: layerName});
        const curlayer = this.ol.getLayerById(layerID);
        if (curlayer)
          this.ol.olCurrentLayer = curlayer;

        //TODO: update currentLayer property when clicked in list (now we hane in issue in updating property)
      };
      this.divLayersList.append(divLayer);
    }
    const layersProperty = this.getProperty('currentLayer');
    if (layersProperty)
      layersProperty.choices = arrLayerNames;
  }

  showCoordsInStatus(p: OLCallbackParam): void {
    if (this.lblStatusCoords) {
      if (p)
        this.lblStatusCoords.innerText = (p.coord[0]).toFixed(3) + ' / ' + (p.coord[1]).toFixed(3);
    }
  }

  showMarkerTooltip(p: OLCallbackParam): void {
    if (!this.showTooltip)
      return;
    if (!p)
      return;
    if (!p.features.length) {
      ui.tooltip.hide;
      this.ol.olMap.render();
      return;
    }
    ui.tooltip.hide;

    let xCrd = this.root.getBoundingClientRect().left;
    let yCrd = 0; //this.root.clientTop;
    if (this.viewerContainer) {
      xCrd += p.pixel[0];// + 150;
      yCrd += p.pixel[1] + 60;
    }
    const markerIdx = p.features[0].get('fieldIndex');
    const tipElement = ui.divV([]);
    for (let cl = 0; cl < this.dataFrame.columns.length; cl ++) {
      const vl = this.dataFrame.columns.byIndex(cl).name + ': ' + this.dataFrame.columns.byIndex(cl).get(markerIdx);
      const subElem = ui.div(vl);
      tipElement.append(subElem);
    }

    ui.tooltip.show(tipElement, xCrd, yCrd);
  }

  handlerOnMarkerSelect(p: OLCallbackParam): void {
    if (!p)
      return;
    if (p.features) {
      let idx = 0;
      this.dataFrame.selection.setAll(false, false);
      for (let i = 0; i < p.features.length; i++) {
        idx = p.features[i].get('fieldIndex');
        //We need to search by element index not by coordinates because coords were tranfrormed while mapping
        if (idx !== undefined)
          this.dataFrame.selection.set(idx, true, false);
      }
      if (idx !== undefined) {
        this.dataFrame.selection.set(idx, true, true);
        this.dataFrame.currentRowIdx = idx; //set focus on the last selected item
      }
    }
  }

  handlerOnMapClick(p: OLCallbackParam): void {
    if (this.lblStatusCoords) {
      if (!p)
        return;
      if (p.features) {
        //perform some actions with clicked features
      }
    }
  }

  private rootOnSizeChanged(args: any): void {
    this.ol.olMap.setSize([this.root.clientWidth, this.root.clientHeight]);
    setTimeout( function(m) {m.updateSize();}, 200, this.ol.olMap);
  }

  onTableAttached(): void {
    this.init();

    if (this.latitudeColumnName === null && this.longitudeColumnName === null) {
      let col = this.dataFrame.columns.bySemType(DG.SEMTYPE.LATITUDE);
      if (col)
        this.latitudeColumnName = col.name;
      col = this.dataFrame.columns.bySemType(SEMTYPEGIS.LATIITUDE);
      if ((col) && (this.latitudeColumnName === null))
        this.latitudeColumnName = col.name;
      col = this.dataFrame.columns.bySemType(DG.SEMTYPE.LONGITUDE);
      if (col)
        this.longitudeColumnName = col.name;
      col = this.dataFrame.columns.bySemType(SEMTYPEGIS.LONGITUDE);
      if ((col) && (this.longitudeColumnName === null))
        this.longitudeColumnName = col.name;
    }

    //events of dataframe handling
    //rows selecting
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 100).subscribe((_) => {
      const selcount = this.dataFrame.selection.getSelectedIndexes();
      if (!this.ol.olMarkersSelLayerGL)
        return;
      if (this.ol.olMarkersSelSource)
        this.ol.olMarkersSelSource.clear();
      for (let i = 0; i < selcount.length; i++) {
        this.ol.addPointSc(this.coordinates[selcount[i]], this.sizeValues[selcount[i]], this.colorValues[selcount[i]],
          this.labelValues[selcount[i]], selcount[i], this.ol.olMarkersSelLayerGL);
      }
      if (!this.ol.preventFocusing) {
        this.ol.olMap.getView().fit(this.ol.olMarkersSelSource.getExtent(),
          {padding: [50, 50, 50, 50], maxZoom: 9});
      }
      this.ol.preventFocusing = false;
      // this.render();
    }));
    //update on filtration
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => {
      this.render(true);
    }));

    this.render(true);
  }

  updateOpenLayerProperties(updateLayer: boolean = false) {
    this.ol.defaultColor = this.defaultColor;
    this.ol.selectedColor = this.selectedColor;
    this.ol.markerMinColor = this.markerMinColor;
    this.ol.markerMaxColor = this.markerMaxColor;
    this.ol.markerOpacity = this.markerOpacity / 100;
    this.ol.markerMinSize = this.markerMinSize;
    this.ol.markerMaxSize = this.markerMaxSize;
    this.ol.markerDefaultSize = this.markerDefaultSize;

    this.ol.heatmapRadius = this.heatmapRadius;
    this.ol.heatmapBlur = this.heatmapBlur;

    if (!this.initialized)
      return;

    if (updateLayer)
      this.ol.updateMarkersGLLayer();

    //TODO: experiment: i've try to refresh WebGL layer here without recreating of it but it's not work yet>>
    // this.ol.olMap.render();
    // window.requestAnimationFrame(() => {
    //   this.ol.olMarkersLayerGL!.updateStyleVariables({'size': this.markerDefaultSize});
    //   this.ol.olMap.render();
    // });
  }

  onPropertyChanged(prop: DG.Property): void {
    if (!this.initialized)
      return;
    if (prop.name === 'currentLayer') {
      this.updateLayersList();
      return;
    }
    if (prop.name === 'heatmapRadius') {
      this.ol.heatmapRadius = this.heatmapRadius;
      return;
    }
    if (prop.name === 'showTooltip')
      return;
    if (prop.name === 'heatmapBlur') {
      this.ol.heatmapBlur = this.heatmapBlur;
      return;
    }
    if ((prop.name === 'defaultColor') ||
    (prop.name === 'selectedColor') ||
    (prop.name === 'markerOpacity') ||
    (prop.name === 'markerMinSize') ||
    (prop.name === 'markerMaxSize') ||
    (prop.name === 'markerDefaultSize')) {
      this.updateOpenLayerProperties(true);
      return;
    }

    this.render();
  }

  detach(): void {
    if (this.ol) {
      this.ol.removeAllLayers();
      this.ol.olMap.dispose();
    }
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  render(fit: boolean = false): void {
    // this.updateOpenLayerProperties(false);

    if (this.latitudeColumnName == null || this.longitudeColumnName == null)
      return;

    this.getCoordinates(); //TODO: not only coordinates but a data at all (?)

    if (this.renderType === 'heat map') {
      this.updateOpenLayerProperties(false);
      this.renderHeat(this.features);
    } else if (this.renderType === 'markers') {
      this.renderMarkersBatch(this.features);
      this.updateOpenLayerProperties(true);
      // this.renderMarkers();
    }

    if (fit) {
      if (this.ol.olMarkersLayerGL)
        this.ol.olMap.getView().fit((this.ol.olMarkersLayerGL.getSource()!).getExtent());
      // this.ol.olMap.getView().fit((this.ol.olMarkersLayer.getSource()!).getExtent());
    }

    this.updateLayersList();
  }

  getCoordinates(): void {
    this.coordinates.length = 0;
    this.labelValues.length = 0;
    this.sizeValues.length = 0;
    this.colorValues.length = 0;
    this.indexValues.length = 0;

    const indexes = this.dataFrame.filter.getSelectedIndexes();
    let labelVal: Int32Array | Float32Array | Float64Array | Uint32Array;
    let colorVal: Int32Array | Float32Array | Float64Array | Uint32Array;
    let sizeVal: Int32Array | Float32Array | Float64Array | Uint32Array;
    //scaling parameters for color and size (TODO: decide wether we need them here or all scaling will be in styles fn)
    let colorCoeff: number;
    let sizeCoeff: number;
    let colorShift: number;
    let sizeShift: number;

    const lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    const lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();
    // val = this.dataFrame.getCol(this.labelsColumnName).getRawData();

    if ((!lat) || (!lon))
      return;

    //TODO: change it to filling array of objects with all table data (if we need it of course?)
    let colValue: DG.Column | null = null;
    let colColor: DG.Column | null = null;
    let colSize: DG.Column | null = null;
    try {
      //TODO: check does it exist corresponding column to prevent exception raising
      if (this.colorColumnName !== '')
        colColor = this.dataFrame.getCol(this.colorColumnName);
      if (this.sizeColumnName !== '')
        colSize = this.dataFrame.getCol(this.sizeColumnName);
      if (this.labelsColumnName !== '')
        colValue = this.dataFrame.getCol(this.labelsColumnName);
    } finally {
      if (colValue)
        labelVal = colValue.getRawData();
      else {
        labelVal = new Float32Array(lat.length); //create array of length equal to latitude array length
        labelVal.fill(1);
      }
      if ((colColor)) {
        //get column color coding settings from column
        this.ol.useColorField = true;
        if (colColor.colors.getType() === DG.COLOR_CODING_TYPE.LINEAR) {
          this.ol.colorCodingType = DG.COLOR_CODING_TYPE.LINEAR;
          this.ol.colorCodingData = colColor.getTag(DG.TAGS.COLOR_CODING_LINEAR);
        } else if (colColor.colors.getType() === DG.COLOR_CODING_TYPE.CATEGORICAL) {
          this.ol.colorCodingType = DG.COLOR_CODING_TYPE.CATEGORICAL;
          this.ol.colorCodingData = colColor.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL);
        } else if (colColor.colors.getType() === DG.COLOR_CODING_TYPE.CONDITIONAL) {
          this.ol.colorCodingType = DG.COLOR_CODING_TYPE.CONDITIONAL;
          this.ol.colorCodingData = colColor.getTag(DG.TAGS.COLOR_CODING_CONDITIONAL);
        } else if ((colColor.type != DG.COLUMN_TYPE.BIG_INT) && (colColor.type != DG.COLUMN_TYPE.INT) &&
                  (colColor.type != DG.COLUMN_TYPE.DATE_TIME) && (colColor.type != DG.COLUMN_TYPE.FLOAT))
          this.ol.useColorField = false;
        // alert(this.ol.colorCodingData);
        //
        colorVal = colColor.getRawData();
        this.ol.minFieldColor = colColor.min;
        this.ol.maxFieldColor = colColor.max;
        if (colColor.max === colColor.min)
          colorCoeff = 1;
        // else colorCoeff = (this.markerMaxColor-this.markerMinColor)/(colColor.max-colColor.min);
        else
          colorCoeff = (0xff0000-0x0000ff) / (colColor.max - colColor.min);
        colorShift = colColor.min;
        //TODO: use another scaling scheme for color (though gradient)
      } else {
        this.ol.useColorField = false;
        colorCoeff = 1;
        colorShift = 0;
        colorVal = new Float32Array(lat.length); //create array of length equal to latitude array length
        colorVal.fill(this.defaultColor);
      }
      //marker size data loading and scale calculation
      //TODO: add checking for number column type
      if ((colSize)) {
        sizeVal = colSize.getRawData();
        this.ol.minFieldSize = colSize.min;
        this.ol.maxFieldSize = colSize.max;
        this.ol.useSizeField = true;
        if (colSize.max === colSize.min)
          sizeCoeff = 1;
        else
          sizeCoeff = Math.abs(this.markerMaxSize - this.markerMinSize) / (colSize.max - colSize.min);
        sizeShift = colSize.min;
      } else {
        this.ol.useSizeField = false;
        sizeCoeff = 1;
        sizeShift = 0;
        sizeVal = new Float32Array(lat.length); //create array of length equal to latitude array length
        sizeVal.fill(this.markerDefaultSize);
      }

      this.features.length = 0; //clear array of features
      for (let i = 0; i < indexes.length; i++) {
        if ((i < lat.length) && (i < lon.length)) {
          this.coordinates.push([lon[indexes[i]], lat[indexes[i]]]);
          this.labelValues.push(labelVal[indexes[i]]);
          this.sizeValues.push(sizeVal[indexes[i]]);
          this.colorValues.push(colorVal[indexes[i]]);
          this.indexValues.push(indexes[i]);
          //this is the old way where we scaled values initially>>
          // this.sizeValues.push((sizeVal[indexes[i]] - sizeShift) * sizeCoeff);
          // this.colorValues.push((colorVal[indexes[i]] - colorShift) * colorCoeff);

          // if (!colColor) this.colorValues.push(this.defaultColor);
          // colColor.meta.markers
          const coords = OLProj.fromLonLat([lon[indexes[i]], lat[indexes[i]]]);
          this.features.push(new Feature({
            geometry: new Point(coords),
            fieldLabel: labelVal[indexes[i]],
            fieldSize: (sizeVal[indexes[i]]),
            fieldColor: (colorVal[indexes[i]]),
            fieldIndex: indexes[i],
            //this is the old way where we scaled values initially>>
            // fieldSize: ((sizeVal[indexes[i]] - sizeShift) * sizeCoeff),
            // fieldColor: toStringColor(((colorVal[indexes[i]] - colorShift) * colorCoeff), this.markerOpacity / 100),
          }));
        }
      }
    } //end of finally block
    //<<getCoordinates function
  }

  renderHeat(arrFeatures: Array<Feature>): void {
    this.getCoordinates();
    let colName = this.colorColumnName;
    if (colName === 'null')
      colName = this.dataFrame.name;
    if (colName === null)
      colName = this.dataFrame.name;

    let layer = this.ol.getLayerByName('HL: ' + colName);
    if (!layer)
      layer = this.ol.addNewHeatMap('HL: ' + colName);

    this.ol.clearLayer(layer);
    this.ol.addFeaturesBulk(arrFeatures, layer);
    // for (let i = 0; i < this.coordinates.length; i++)
    //   this.ol.addPointSc(this.coordinates[i], this.sizeValues[i], this.colorValues[i], this.values[i], layer);
  }

  renderMarkers(): void {
    this.ol.clearLayer(); //TODO: clear exact layer
    for (let i = 0; i < this.coordinates.length; i++) {
      this.ol.addPointSc(this.coordinates[i], this.sizeValues[i], this.colorValues[i],
        this.labelValues[i], i);
    }
  }

  renderMarkersBatch(arrFeatures: Array<Feature>): void {
    this.ol.clearLayer(); //TODO: clear exact layer
    this.ol.addFeaturesBulk(arrFeatures);
  }
}
