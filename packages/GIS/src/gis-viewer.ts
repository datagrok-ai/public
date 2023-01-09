/*
 * GIS VIEWER
 * this is the viewer module which manage viewer UI, map based on open-layers class
*/
//Base import
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//GIS semantic types import
import {SEMTYPEGIS} from '../src/gis-semtypes';

//OpenLayers functionality import
import {OLCallbackParam, OpenLayers, Coordinate, toStringColor} from '../src/gis-openlayer';
import VectorLayer from 'ol/layer/Vector';
import * as OLProj from 'ol/proj';
import Feature from 'ol/Feature';
//import {Circle} from 'ol/geom';
import {Point} from 'ol/geom';
// import { numberSafeCompareFunction } from 'ol/array';

//release mode flag: set to true to hide all experimental UI features
const releaseMode = false;

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
  renderType: string;
  heatmapRadius: number;
  heatmapBlur: number;
  showTooltip: boolean;
  autoScaling: boolean;

  isShortUI: boolean = true;

  //ui elements
  panelLeft: HTMLElement | null = null;
  btnLeftPanel: HTMLElement | null = null;
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
  // layers = [];
  features: Array<Feature> = [];
  featuresFull: Array<Feature> = [];
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
    this.currentLayer = this.string('currentLayer', 'BaseLayer', {category: 'Layers', userEditable: false});
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    this.colorColumnName = this.string('colorColumnName');
    this.sizeColumnName = this.string('sizeColumnName');
    this.labelsColumnName = this.string('labelsColumnName', '', {userEditable: false});
    this.markersColumnName = this.string('markersColumnName', '', {userEditable: false});

    this.markerOpacity = this.float('markerOpacity', 80, {category: 'Markers', editor: 'slider', min: 0, max: 100});
    this.markerDefaultSize = this.int('markerDefaultSize', 5, {category: 'Markers', min: 1, max: 30});
    //DONE: there is no need in valitators: Min can be > Max
    this.markerMinSize = this.int('markerMinSize', 2, {category: 'Markers', min: 1, max: 30});
    this.markerMaxSize = this.int('markerMaxSize', 15, {category: 'Markers', min: 1, max: 30});

    this.defaultColor = this.int('defaultColor', 0x1f77b4, {category: 'Markers'});
    this.selectedColor = this.int('selectedColor', 0xff8c00, {category: 'Markers'});
    this.markerMinColor = this.int('markerMinColor', 0x0000ff, {category: 'Markers', userEditable: false});
    this.markerMaxColor = this.int('markerMaxColor', 0xff0000, {category: 'Markers', userEditable: false}); //hexToRGB

    this.heatmapRadius = this.int('heatmapRadius', 3,
      {category: 'Heatmap', description: 'Heatmap radius', min: 1, max: 10, userEditable: true});
    this.heatmapBlur = this.int('heatmapBlur', 30,
      {category: 'Heatmap', description: 'Heatmap radius', min: 1, max: 50, userEditable: true});

    this.showTooltip = this.bool('showTooltip', false);
    this.autoScaling = this.bool('autoScaling', false);
    this.renderType = this.string('renderType', 'markers', {choices: ['markers', 'heatmap', 'both']});
    //<<end of constructor()
  }

  layerUIElement(num: number, layerName: string, layerId: string,
    isVisible: boolean, isExportable: boolean = true): HTMLElement {
    //
    let btnExport: HTMLElement = ui.div();
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
    //Visible button>>
    const btnVisible = ui.iconFA(isVisible ? 'eye' : 'eye-slash', (evt) => {
      evt.stopImmediatePropagation();
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
    setupBtn(btnVisible, layerName, layerId);

    //Export button>>
    if (isExportable) {
      btnExport = ui.iconFA('download', (evt) => { //cogs
        evt.stopImmediatePropagation();
        const divLayer = (evt.currentTarget as HTMLElement);
        const layerId = divLayer.getAttribute('layerId');
        if (!layerId)
          return;
        const layer = this.ol.getLayerById(layerId) as VectorLayer<any>;
        if (!layer)
          return;
        const arrPreparedToDF: any[] = this.ol.exportLayerToArray(layer);
        if (arrPreparedToDF.length) {
          const df = DG.DataFrame.fromObjects(arrPreparedToDF);
          if (df) {
            const gisCol = df.col('gisObject');
            if (gisCol)
              gisCol.semType = SEMTYPEGIS.GISAREA; //SEMTYPEGIS.GISOBJECT;
            df.name = layer.get('layerName');

            // this.view.addTableView(df as DG.DataFrame);
            grok.shell.addTableView(df as DG.DataFrame);
          }
        }
      }, 'Export layer data to table');
      setupBtn(btnExport, layerName, layerId);
    } //<<if the flag isExportable

    //Delete button>>
    const btnDelete = ui.iconFA('trash-alt', (evt)=>{ //trash
      const divLayer = (evt.currentTarget as HTMLElement);
      const layerId = divLayer.getAttribute('layerId');
      if (layerId)
        this.ol.removeLayerById(layerId);
      this.updateLayersList();
    }, 'Delete layer');
    setupBtn(btnDelete, layerName, layerId);

    const panelButtons = ui.divH([btnVisible, btnDelete]);

    //add export button if the flag isExportable == true
    if (isExportable)
      panelButtons.append(btnExport);

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
    if ((this.panelLeft) && (this.btnLeftPanel)) {
      this.panelLeft.style.visibility = shortUI ? 'visible' : 'hidden';
      this.btnLeftPanel.click();
    }
  }

  initUi(shortUI: boolean = true): HTMLElement {
    this.isShortUI = shortUI;
    //create UI>>
    const uiButtonsToAdd: HTMLElement[] = [];

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
    this.btnLeftPanel = ui.button(ui.iconFA('layer-group'), ()=>{ //cogs
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
    }, 'Show/hide layers panel');
    this.btnLeftPanel.style.margin = '1px';
    uiButtonsToAdd.push(this.btnLeftPanel);

    //Heatmap button
    const btnHeatmap = ui.button('HM', () => {
      this.getCoordinates();
      this.renderHeat(this.features);
      this.updateLayersList();
    }, 'Build heatmap for data');
    btnHeatmap.style.margin = '1px';
    uiButtonsToAdd.push(btnHeatmap);

    //Export selected button
    const btnExportSelected = ui.button(ui.iconFA('copy'), () => {
      // this.dataFrame.selection
      const selDF = this.dataFrame.clone(this.dataFrame.selection);
      selDF.name = this.dataFrame.name + ' selection';
      grok.shell.addTableView(selDF as DG.DataFrame);
    }, 'Export selected markers to new table');
    btnExportSelected.style.margin = '1px';
    if (!releaseMode)
      uiButtonsToAdd.push(btnExportSelected);

    //Select features by area button
    const btnSelectByArea = ui.button(ui.iconFA('lasso'), () => {
      if (!this.ol)
        return;
      if (!this.ol.currentAreaObject)
        return;
      this.ol.selectMarkersByGeometry(this.ol.currentAreaObject.getGeometry());
    }, 'Select markers by area');
    btnSelectByArea.style.margin = '1px';
    if (!releaseMode)
      uiButtonsToAdd.push(btnSelectByArea);

    //Select features by area button
    const btnSelectByDistance = ui.button(ui.iconFA('bullseye'), () => {
      //temp testing code>>
      this.ol.testFunc();

      // if (!this.ol)
      //   return;
      // if (!this.ol.currentAreaObject)
      //   return;
      // if (this.ol.currentAreaObject.getGeometry()?.getType() === 'Point') {
      //   const coords = (this.ol.currentAreaObject.getGeometry() as Point).getCoordinates();
      //   const searchArea = new Circle(coords, 500);
      //   this.ol.selectMarkersByGeometry(searchArea);
      // }

      // const searchFT = new Feature({
      //   geometry: new Circle(coords, ) };
    }, 'Select within distance');
    btnSelectByDistance.style.margin = '1px';
    if (!releaseMode)
      uiButtonsToAdd.push(btnSelectByDistance);

    //add buttons to top menu panel
    this.panelTop.append(ui.divH(uiButtonsToAdd));
    // this.panelTop.append(ui.divH([this.btnLeftPanel, btnHeatmap, btnExportSelected]));
    // this.panelTop.append(ui.divH([this.btnLeftPanel, heatmapBtn]));

    //left panel icons>>
    this.divLayersList = ui.divV([]);
    this.panelLeft.append(this.divLayersList);

    this.panelLeft.style.visibility = 'visible';
    this.btnLeftPanel.click();

    const body = ui.box();
    body.id = 'map-container';
    body.style.maxWidth = '100%';
    body.style.maxHeight = '100%';

    this.viewerContainer = body;/*ui.splitV(
      [this.panelTop,
        ui.splitH([this.panelLeft, body], null, true),
        this.panelBottom]);*/

    this.root.appendChild(this.viewerContainer);

    //setup context menu
    this.onContextMenu.subscribe((menu) => {
      // menu.item('Add rest layer', async () => {
      //   this.ol.testFunc();
      // });
      menu.item('Reset view', async () => {
        await this.render(true, false);
      });
    });

    return this.viewerContainer;
  }

  init() {
    const loadingDiv = ui.div('Loading...');
    try {
      loadingDiv.style.position = 'absolute';
      loadingDiv.style.left = '50%';
      loadingDiv.style.top = '50%';

      this.root.append(loadingDiv);

      this.ol.useWebGL = true; //choose true if we want to use WebGL renderer
      this.initUi();
      this.ol.initMap('map-container');

      this.updateOpenLayerProperties();

      this.updateLayersList();

      //subscribe to events
      this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
      this.subs.push(ui.onSizeChanged((this.panelLeft as HTMLElement)).subscribe(this.rootOnSizeChanged.bind(this)));
      //setup callbacks
      this.ol.setMapPointermoveCallback(this.showMarkerTooltip.bind(this));
      this.ol.setMapSelectionCallback(this.handlerOnMarkerSelect.bind(this));
      this.ol.setMapClickCallback(this.handlerOnMapClick.bind(this));
      this.ol.setMapRefreshCallback(this.updateLayersList.bind(this));

      this.initialized = true;
    } catch (e: any) {
      this.initialized = false;
      grok.shell.error(e.toString());
    } finally {
      setTimeout(() => {
        // const loadingDiv = this.root.firstElementChild as HTMLElement;
        //loadingDiv.style.visibility = 'hidden';
        this.root.removeChild(loadingDiv);

        grok.shell.o = this;
        grok.shell.windows.showProperties = true;
      }, 100);
    }
  }

  updateLayersList(): void {
    if (this.isShortUI)
      return;
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

        //TODO: update currentLayer property when clicked in list (now we have an issue in updating properties)
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
    if ((!this.showTooltip) || (!p))
      return;

    ui.tooltip.hide;
    this.ol.olMap.render();
    if (!p.features.length) {
      // ui.tooltip.hide;
      // this.ol.olMap.render();
      return;
    }

    let xCrd = this.root.getBoundingClientRect().left;
    let yCrd = 0;
    if (this.viewerContainer) {
      xCrd += p.pixel[0];
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
        //We need to search by element index because coords were transformed while mapping
        if (idx !== undefined)
          this.dataFrame.selection.set(idx, true, false);
      }
      this.dataFrame.selection.fireChanged();
      if (idx !== undefined)
        this.dataFrame.currentRowIdx = idx; //set focus on the last selected item
    }
  }

  handlerOnMapClick(p: OLCallbackParam): void {
    if (this.lblStatusCoords) {
      if (!p)
        return;
      if (p.features && p.features.length > 0) {
        //perform some actions with clicked features
      } else { //if no features were clicked - just show properties for map
        setTimeout(() => {
          grok.shell.o = this;
          grok.shell.windows.showProperties = true;
        }, 50);
      }
    }
  }

  private rootOnSizeChanged(args: any): void {
    this.ol.olMap.setSize([this.root.clientWidth, this.root.clientHeight]);
    setTimeout( function(m) {m.updateSize();}, 200, this.ol.olMap);
  }

  async onTableAttached(): Promise<void> {
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
        this.ol.addPoint(this.coordinates[selcount[i]], this.sizeValues[selcount[i]], this.colorValues[selcount[i]],
          this.labelValues[selcount[i]], selcount[i], this.ol.olMarkersSelLayerGL);
      }
      if ((!this.ol.preventFocusing) && (selcount.length > 0)) {
        this.ol.olMap.getView().fit(this.ol.olMarkersSelSource.getExtent(),
          {padding: [50, 50, 50, 50], maxZoom: 9});
      }
      this.ol.preventFocusing = false;
    }));

    //update on filtration
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 200).subscribe((_) => {
      const indexes = this.dataFrame.filter.getSelectedIndexes();
      //new way of filtration
      this.features.length = 0;
      for (let i = 0; i < indexes.length; i++)
        this.features.push(this.featuresFull[indexes[i]]);

      this.render(this.autoScaling, false);

      //old way of filtration>>
      //Unfortunately it doesn't work for the heatmap layer
      // for (let i = 0; i < this.features.length; i++)
      //   this.features[i].set('filtered', 0);
      // for (let i = 0; i < indexes.length; i++)
      //   this.features[indexes[i]].set('filtered', 1);
      //<<old way of filtration

      // this.ol.olMarkersSource.changed();
      // this.ol.updateMarkersGLLayer(false);
      // this.ol.olMap.render();
    }));

    //update on color coding change
    this.subs.push(DG.debounce(this.dataFrame.onMetadataChanged, 100).subscribe((data) => {
      const eventColumn = ((data.args.source as unknown) as DG.Column);
      if (!eventColumn)
        return;

      // eventColumn.meta.colors.getColors();
      if (eventColumn.name === this.colorColumnName) {
        this.refreshColorCodingStyle(eventColumn); //refilling of color coding
        this.ol.updateMarkersGLLayer();
        // this.ol.olMarkersLayerGL?.changed();
        // this.ol.olMarkersSource.changed();
        // this.ol.olMap.render();
      }
    }));

    //render map
    setTimeout(async () => {
      await this.render(true, true);
    }, 2);

  //<< end of onTableAttached method
  }

  updateOpenLayerProperties(updateLayer: boolean = false): void {
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

    if (updateLayer) {
    //TODO: experiment: i've tried to refresh WebGL layer here without recreating of it but it isn't work yet>>
      this.ol.updateMarkersGLLayer(true);
      // this.ol.olMarkersSource.changed();
      // this.ol.olMap.render();
    }
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
    if (prop.name === 'heatmapBlur') {
      this.ol.heatmapBlur = this.heatmapBlur;
      return;
    }
    if (prop.name === 'showTooltip')
      return;
    if (prop.name === 'autoScaling')
      return;

    if ((prop.name === 'defaultColor') ||
    (prop.name === 'selectedColor') ||
    (prop.name === 'markerOpacity') ||
    (prop.name === 'markerMinSize') ||
    (prop.name === 'markerMaxSize') ||
    (prop.name === 'markerDefaultSize')) {
      this.updateOpenLayerProperties(true);
      return;
    }

    this.render(this.autoScaling, true);
  }

  detach(): void {
    if (this.ol) {
      this.ol.removeAllLayers();
      this.ol.olMap.dispose();
    }
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async render(fit: boolean = false, reloadData: boolean = true): Promise<void> {
    try {
      if (reloadData)
        this.getCoordinates();

      this.updateOpenLayerProperties(true);
      if (this.renderType === 'heatmap') {
        //render Heatmap
        this.renderHeat(this.features);
        //TODO: this style of switching visibility is a bad but temporary decision
        this.ol.olHeatmapLayer?.setVisible(true);
        this.ol.olMarkersLayerGL?.setVisible(false);
      } else if (this.renderType === 'markers') {
        //render markers map
        this.renderMarkersBatch(this.features);
        this.ol.olHeatmapLayer?.setVisible(false);
        this.ol.olMarkersLayerGL?.setVisible(true);
      } else if (this.renderType === 'both') {
        //render markers map
        this.renderMarkersBatch(this.features);
        this.ol.olHeatmapLayer?.setVisible(true);
        this.ol.olMarkersLayerGL?.setVisible(true);
      }

      if (fit) {
        if ((this.ol.olMarkersLayerGL) && (this.features.length > 0))
          this.ol.olMap.getView().fit((this.ol.olMarkersSource).getExtent());
      }
    } finally {
      this.updateLayersList();
    }
  }

  refreshColorCodingStyle(colColor: DG.Column): void {
    this.ol.useColorField = true;
    this.ol.colorCodingData = '';
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
    let colorCodes: Uint32Array = new Uint32Array();

    if (!this.dataFrame)
      return;

    let lat = null;
    let lon = null;
    lat = this.dataFrame.col(this.latitudeColumnName)?.getRawData();
    lon = this.dataFrame.col(this.longitudeColumnName)?.getRawData();

    // clear the features full array
    if ((!lat) || (!lon)) {
      this.featuresFull.length = 0;
      this.features.length = 0;
      return;
    }


    //TODO: change it to filling array of objects with all table data (if we need it of course?)
    let colValue: DG.Column | null = null;
    let colColor: DG.Column | null = null;
    let colSize: DG.Column | null = null;
    try {
      //DONE: check does it exist corresponding column to prevent exception raising
      if (this.colorColumnName !== '')
        colColor = this.dataFrame.col(this.colorColumnName);
      if (this.sizeColumnName !== '')
        colSize = this.dataFrame.col(this.sizeColumnName);
      if (this.labelsColumnName !== '')
        colValue = this.dataFrame.col(this.labelsColumnName);
    } finally {
      if (colValue)
        labelVal = colValue.getRawData(); //!!!
      else {
        labelVal = new Float32Array(lat.length); //create array of length equal to latitude array length
        labelVal.fill(0);
      }
      if ((colColor)) {
        //get column color coding settings from column
        this.refreshColorCodingStyle(colColor);

        //new way of color coding
        colorCodes = colColor.meta.colors.getColors();
        if (colorCodes.length === 0) {
          colorVal = new Uint32Array(lat.length);
          colorCodes.fill(this.defaultColor);
        }

        colorVal = colColor.getRawData();
        this.ol.minFieldColor = colColor.min;
        this.ol.maxFieldColor = colColor.max;
      } else {
        this.ol.useColorField = false;
        colorVal = new Float32Array(lat.length); //create array of length equal to latitude array length
        colorVal.fill(this.defaultColor);
      }
      //marker size data loading and scale calculation
      if ((colSize) &&
      ((colSize.type === DG.COLUMN_TYPE.BIG_INT) || (colSize.type === DG.COLUMN_TYPE.INT) ||
      (colSize.type === DG.COLUMN_TYPE.DATE_TIME) || (colSize.type === DG.COLUMN_TYPE.FLOAT)) ) {
        sizeVal = colSize.getRawData();
        this.ol.minFieldSize = colSize.min;
        this.ol.maxFieldSize = colSize.max;
        this.ol.useSizeField = true;
      } else {
        this.ol.useSizeField = false;
        sizeVal = new Float32Array(lat.length); //create array of length equal to latitude array length
        sizeVal.fill(this.markerDefaultSize);
      }

      this.featuresFull.length = 0;
      for (let i = 0; i < lat.length; i++) {
        // if ((lon[i] !== DG.FLOAT_NULL) && (lat[i] !== DG.FLOAT_NULL))
        {
          this.coordinates.push([lon[i], lat[i]]);
          this.labelValues.push(labelVal[i]);
          this.sizeValues.push(sizeVal[i]);
          this.colorValues.push(colorVal[i]);
          this.indexValues.push(i);

          const coords = OLProj.fromLonLat([lon[i], lat[i]]);
          const ft = new Feature({
            geometry: new Point(coords),
            fieldLabel: labelVal[i],
            fieldSize: (sizeVal[i]),
            fieldColor: (colorVal[i]),
            fieldColorCode: (toStringColor(colorCodes[i], this.markerOpacity)),
            fieldIndex: i,
            //TODO: test this approach of filtering-out data with blank fields (for cases of heat-map, etc.)
            filtered: ((lon[i] !== DG.FLOAT_NULL) && (lat[i] !== DG.FLOAT_NULL)) ? 1 : 0,
          });
          this.featuresFull.push(ft);
        }
      }

      this.features.length = 0; //clear array of features
      for (let i = 0; i < indexes.length; i++)
        this.features.push(this.featuresFull[indexes[i]]);
    } //end of finally block
    //<<getCoordinates function
  }

  renderHeat(arrFeatures: Array<Feature>): void {
    // this.getCoordinates();
    let colName = this.colorColumnName;
    if (colName === '')
      colName = this.dataFrame.name;

    //for now we we don't need multiple heatmap layers - so "'HL: ' + colName" code commented
    let layer = this.ol.getLayerByName('Heatmap'); //'HL: ' + colName);
    if (!layer)
      layer = this.ol.addNewHeatMap('Heatmap'); //'HL: ' + colName);

    this.ol.clearLayer(layer);
    this.ol.addFeaturesBulk(arrFeatures, layer);
  }

  renderMarkers(): void {
    this.ol.clearLayer(); //NOTE: for now we have common source for layers (so don't need to search for layer)
    for (let i = 0; i < this.coordinates.length; i++) {
      this.ol.addPoint(this.coordinates[i], this.sizeValues[i], this.colorValues[i],
        this.labelValues[i], i);
    }
  }

  renderMarkersBatch(arrFeatures: Array<Feature>): void {
    this.ol.clearLayer(); //NOTE: for now we have common source for layers (so don't need to search for layer)
    this.ol.addFeaturesBulk(arrFeatures);
  }
}
