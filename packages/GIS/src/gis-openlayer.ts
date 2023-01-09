//base import
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
//import * as ui from 'datagrok-api/ui';

import * as GisTypes from '../src/gis-semtypes';
import {GisViewer} from './gis-viewer';

import {Map as OLMap, MapBrowserEvent, View as OLView} from 'ol';
import HeatmapLayer from 'ol/layer/Heatmap';
import BaseLayer from 'ol/layer/Base';
import Layer from 'ol/layer/Layer';
import TileLayer from 'ol/layer/Tile';
import VectorLayer from 'ol/layer/Vector';
import WebGLPointsLayer from 'ol/layer/WebGLPoints';
// import TileImage from 'ol/source/TileImage'; //this is the base class for XYZ, BingMaps etc..
import LayerRenderer from 'ol/renderer/Layer';
//Sources import
import OSM from 'ol/source/OSM';
import BingMaps from 'ol/source/BingMaps';
import VectorSource from 'ol/source/Vector';
//import XYZ from 'ol/source/XYZ';
import {StyleLike} from 'ol/style/Style';
//Projections working itilities
import * as OLProj from 'ol/proj';
import {Coordinate} from 'ol/coordinate';
//geometry drawing funtions
import Feature, {FeatureLike} from 'ol/Feature';
import Collection from 'ol/Collection';
import * as OLGeom from 'ol/geom';
import * as OLStyle from 'ol/style';
import {LiteralStyle} from 'ol/style/literal';
//import interactions and events
import * as OLInteractions from 'ol/interaction';
import * as OLEventsCondition from 'ol/events/condition';
// import * as OLEvents from 'ol/events';

//import processors
import {GPX, GeoJSON, IGC, KML, TopoJSON} from 'ol/format';
import Source from 'ol/source/Source';
//import {Control} from 'ol/control';
import {defaults as defaultControls} from 'ol/control';

import {PanelLayersControl, BtnLayersControl} from './gis-mapcontrols';

//for test func
import TileWMS from 'ol/source/TileWMS';
//import ImageLayer from 'ol/layer/Image';

//ZIP utilities
import JSZip from 'jszip';
import TileImage from 'ol/source/TileImage';

export {Coordinate} from 'ol/coordinate';

type WebGLPts = WebGLPointsLayer<VectorSource<OLGeom.Point>>;

//interface for callback functions parameter
export interface OLCallbackParam {
  coord: Coordinate; //[number, number];
  pixel: [number, number];
  features: FeatureLike[];
}

let OLG: OpenLayers; //TODO: remove this terrible stuff!
let renderTime: number = 0; //temporary code for benchmarking


export function toStringColor(num : number, opacity?: number) : string {
  num >>>= 0;
  const b = num & 0xFF;
  const g = (num & 0xFF00) >>> 8;
  const r = (num & 0xFF0000) >>> 16;
  const a = opacity ? opacity : 1;
  return 'rgba(' + [r, g, b, a].join(',') + ')';
}

function rgbToHex(color: string): string {
  const rgb = color.split(',');
  const hex = parseInt(rgb[0]).toString(16) + parseInt(rgb[1]).toString(16) + parseInt(rgb[2]).toString(16);
  return hex;
}

// Define a KMZ format class by subclassing ol/format/KML
export async function getKMZData(buffer: any): Promise<string> {
  const zip = new JSZip();
  let kmlData = '';
  await zip.loadAsync(buffer);
  const kmlFile = zip.file(/.kml$/i)[0];
  if (kmlFile)
    kmlData = await kmlFile.async('string');
  return kmlData;
}

export class KMZ extends KML {
  constructor(opt: any) {
    const options = opt || {};
    // options.iconUrlFunction = getKMLImage;
    super(options);
  }
  getType(): any {
    return 'arraybuffer'; // @typedef {'arraybuffer' | 'json' | 'text' | 'xml'} Type
  }
  readFeature(source: any, options: any) {
    const kmlData = getKMZData(source);
    return super.readFeature(kmlData, options);
  }
  readFeatures(source: any, options: any) {
    const kmlData = getKMZData(source);
    return super.readFeatures(kmlData, options);
  }
}

//functions to benchmark layer rendering and change render style to improve performance
function onPostrenderMarkers() {
  renderTime = Date.now() - renderTime;
  console.log('Markers render time = ' + renderTime);
}

function onPrerenderMarkers() {
  renderTime = Date.now();
}


export class OpenLayers {
  olMap: OLMap;
  olCurrentView: OLView;
  olCurrentLayer: BaseLayer | VectorLayer<VectorSource> | HeatmapLayer | WebGLPts | null;
  olBaseLayer: BaseLayer | null;
  olMarkersLayer: VectorLayer<VectorSource> | null;
  olMarkersLayerGL: WebGLPts | null;
  olMarkersSelLayerGL: WebGLPts | null;
  olHeatmapLayer: HeatmapLayer | null;
  //vector sources for markers (now we have just one source for Markers and Heatmap layers)
  olMarkersSource: VectorSource<OLGeom.Point>;
  olMarkersSelSource: VectorSource<OLGeom.Point>;
  olSelectedMarkers: Collection<Feature>;

  currentAreaObject: Feature | null = null;

  useWebGLFlag: boolean = true;
  preventFocusing: boolean = false;
  //map interacions
  dragAndDropInteraction: OLInteractions.DragAndDrop;
  selectInteraction: OLInteractions.Select;
  dragBox: OLInteractions.DragBox;

  //map controls
  btnLayersControl: BtnLayersControl;
  panelLayersList: PanelLayersControl;

  //map event handlers
  onSelectCallback: Function | null = null;
  onClickCallback: Function | null = null;
  onPointermoveCallback: Function | null = null;
  onRefreshCallback: Function | null = null;

  //styles>>
  styleVectorLayer: OLStyle.Style;
  styleVectorSelLayer: OLStyle.Style;

  markerGLStyle: LiteralStyle;
  markerGLSelStyle: LiteralStyle;

  //properties from viewer
  markerDefaultSize: number = 3;
  markerMinSize: number = 2;
  markerMaxSize: number = 15;
  markerOpacity: number = 0.8;
  defaultColor: number = 0x1f77b4;
  selectedColor: number = 0xff8c00;
  markerMinColor: number = 0x0000ff;
  markerMaxColor: number = 0xff0000;
  heatmapBlurParam: number = 20;
  heatmapRadiusParam: number = 10;

  //boundary properties for fieldSize, fieldColor (to interpolate in WebGL style function)
  minFieldColor: number = 0;
  maxFieldColor: number = 100;
  useColorField: boolean = false;
  minFieldSize: number = 0;
  maxFieldSize: number = 100;
  useSizeField: boolean = false;
  colorCodingType: DG.ColorCodingType = DG.COLOR_CODING_TYPE.OFF;
  colorCodingData: string = '';

  constructor(gV?: GisViewer) {
    this.olMap = new OLMap({});
    this.olCurrentView = new OLView({});
    this.olCurrentLayer = null;
    this.olBaseLayer = null;
    this.olMarkersLayer = null;
    this.olMarkersLayerGL = null;
    this.olMarkersSelLayerGL = null;
    this.olHeatmapLayer = null;
    this.markerGLStyle = {};
    this.markerGLSelStyle = {};
    this.olSelectedMarkers = new Collection<Feature>;

    this.olMarkersSource = new VectorSource<OLGeom.Point>();
    this.olMarkersSelSource = new VectorSource<OLGeom.Point>();

    //controls
    this.btnLayersControl = new BtnLayersControl(null);
    this.panelLayersList = new PanelLayersControl(this, null);

    this.styleVectorLayer = new OLStyle.Style({
      fill: new OLStyle.Fill({
        color: '#eeeeee',
      }),
      stroke: new OLStyle.Stroke({
        color: 'rgba(200, 0, 0, 1)',
        width: 1,
      }),
    });
    this.styleVectorSelLayer = new OLStyle.Style({
      fill: new OLStyle.Fill({
        color: 'rgba(200, 50, 50, 0.5)',
      }),
      stroke: new OLStyle.Stroke({
        color: 'rgba(255, 0, 0, 1)',
        width: 3,
      }),
    });

    this.dragAndDropInteraction = new OLInteractions.DragAndDrop({
      formatConstructors: [
        KML,
        // KMZ,
        GPX,
        GeoJSON,
        IGC,
        TopoJSON,
      ],
    });
    this.dragAndDropInteraction.on('addfeatures', this.dragNdropInteractionFn.bind(this));

    this.selectInteraction = new OLInteractions.Select({
      condition: OLEventsCondition.click,
      multi: true,
      style: this.styleVectorSelLayer,
      layers: this.selectCondition.bind(this),
      features: this.olSelectedMarkers,
    });

    this.selectInteraction.on('select', (e) => {
      this.updateSelection(e.target.getFeatures());
    });

    this.dragBox = new OLInteractions.DragBox({
      condition: OLEventsCondition.platformModifierKeyOnly,
    });

    this.dragBox.on('boxstart', (e) => {
      this.selectInteraction.getFeatures().clear();
    });

    this.dragBox.on('boxend', () => {
      this.selectMarkersByGeometry(this.dragBox.getGeometry());
    });

    OLG = this;
    //<<end of constructor
  }

  set useWebGL(v: boolean) {
    this.useWebGLFlag = v;
    this.olMarkersLayer?.setVisible(!this.useWebGLFlag);
    this.olMarkersLayerGL?.setVisible(this.useWebGLFlag);
    this.olMarkersSelLayerGL?.setVisible(this.useWebGLFlag);
  }
  get useWebGL(): boolean {return this.useWebGLFlag;}

  set heatmapBlur(val: number) {
    this.heatmapBlurParam = val;
    if (this.olHeatmapLayer)
      this.olHeatmapLayer.setBlur(this.heatmapBlurParam);
    //TODO: search and apply parameter to the first heatmap layer (if this is needed)
    //NOTE: this will be usefull in case of multilayers approach - for a now we have just one layer
  }
  get heatmapBlur(): number {return this.heatmapBlurParam;}

  set heatmapRadius(val: number) {
    this.heatmapRadiusParam = val;
    if (this.olHeatmapLayer)
      this.olHeatmapLayer.setRadius(this.heatmapRadiusParam);
  }
  get heatmapRadius(): number {return this.heatmapRadiusParam;}

  selectCondition(lr: Layer<Source, LayerRenderer<any>>): boolean {
    //here we can detect whether we want to select object on the map
    //for example (commented below) we can choose just objects from layer: olMarkersLayerGL
    return true; //((lr === this.olMarkersLayer) || (lr === this.olMarkersLayerGL));
  }

  testFunc(): void {
    const Tl = new TileLayer({
      source: new TileImage({ //new XYZ({
        // url: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}',
        url: 'https://tile.openweathermap.org/map/wind_new/{z}/{x}/{y}.png?appid=809f4d16303ce4ae52da96eea0fadc6a',
        // url: 'https://sampleserver1.arcgisonline.com/ArcGIS/rest/services/Specialty/ESRI_StateCityHighway_USA/MapServer',
        // url: 'https://ndmc-001.unl.edu:8080/cgi-bin/mapserv.exe?map=/ms4w/apps/usdm/service/usdm_20221213_wms.map&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&LAYERS=usdm20221213&WIDTH=640&HEIGHT=480&crs=EPSG:3857&styles=default&format=image/png&bbox=-18367715.9809,1689200.13961,-6679169.4476,15538711.0963'
      }),
    });
    Tl.set('layerName', 'TestLayer');
    this.addLayer(Tl);

    const Tl2 = new TileLayer({
      // source: new TileWMS({
      //   url: 'https://ahocevar.com/geoserver/wms',
      //   params: {'LAYERS': 'ne:ne'},
      //   serverType: 'geoserver',
      //   crossOrigin: 'anonymous',
      // }),

      // source: new TileWMS({
      //   url: 'https://mesonet.agron.iastate.edu/cgi-bin/wms/nexrad/n0r-t.cgi',
      //   params: {'LAYERS': 'nexrad-n0r-wmst'},
      //   // params: {'LAYERS': 'uscounties'},
      //   serverType: 'geoserver',
      //   crossOrigin: 'anonymous',

      source: new TileWMS({
        url: 'https://ahocevar.com/geoserver/wms',
        params: {'LAYERS': 'topp:states'}, // , 'TILED': true},
        serverType: 'geoserver',
        // Countries have transparency, so do not fade tiles:
        transition: 0,
      }),
    });
    Tl2.set('layerName', 'WMSTestLayer');
    this.addLayer(Tl2);

    return;
  }

  initMap(targetName: string) {
    if (targetName === '')
      return;

    this.olMap = new OLMap({
      target: targetName,
      controls: defaultControls({attribution: false, rotate: false}).extend([
        this.btnLayersControl,
        this.panelLayersList]),
      view: new OLView({
        projection: 'EPSG:3857', // projection: 'EPSG:4326',
        center: OLProj.fromLonLat([-77.01072, 38.91333]),
        zoom: 5,
        enableRotation: false,
      }),
    });

    //add layers>>
    this.addNewBingLayer('Bing sat'); //TODO: if we want to add sattelite layer we have to manage it properly
    this.olBaseLayer = this.addNewOSMLayer('BaseLayer');

    //prepare markersGLLayer
    this.updateMarkersGLLayer();

    //prepare markersLayer
    if (!this.useWebGL) {
      // eslint-disable-next-line max-len
      this.olMarkersLayer = this.addNewVectorLayer('Markers', null, this.genStyleMarker.bind(this), this.olMarkersSource);
      this.olMarkersLayer.on('postrender', onPostrenderMarkers);
      this.olMarkersLayer.on('prerender', onPrerenderMarkers);
    }

    this.useWebGL = true;

    this.olHeatmapLayer = this.addNewHeatMap('Heatmap');
    this.olHeatmapLayer.setVisible(false);

    //add base event handlers>>
    this.olMap.on('click', this.onMapClick.bind(this));
    this.olMap.on('pointermove', this.onMapPointermove.bind(this));

    //add interactions to the map
    this.olMap.addInteraction(this.selectInteraction); //select interaction
    this.olMap.addInteraction(this.dragBox); //add box selection interaction
    this.olMap.addInteraction(this.dragAndDropInteraction); //add dragNdrop ability

    this.setBtnLayersClickCallback(this.onLayersListBtnClick.bind(this));
    this.panelLayersList.setVisibility(false);
    //this.onLayersListBtnClick();
    //<<end of InitMap function
  }

  updateSelection(sel: Collection<Feature<OLGeom.Geometry>>) {
    const res: OLCallbackParam = {
      coord: [0],
      pixel: [0, 0],
      features: [],
    };
    // const src = this.olMarkersSelLayerGL.getSource();
    // sel.forEach((ft: Feature) => {
    //   const geom = ft.getGeometry();
    //   if (geom?.getType() === 'Point')
    //     res.features.push(ft);
    // }); //<<for each feature
    //<<the way commented above is too slow
    //(but in case of pushing the whole array there is a risk to push not only Points)
    res.features.push(...sel.getArray());
    if (this.olMarkersSelLayerGL) {
      this.clearLayer(this.olMarkersSelLayerGL);
      this.addFeaturesBulk((res.features as Feature<OLGeom.Geometry>[]), this.olMarkersSelLayerGL);
      this.olMarkersSelLayerGL.changed();
    }
    this.preventFocusing = true;
    if (this.onSelectCallback)
      this.onSelectCallback(res);
  }

  selectMarkersByGeometry(geom: OLGeom.Geometry | undefined | null): void {
  // selectMarkersByGeometry(geom: OLGeom.Geometry): Collection<Feature<OLGeom.Geometry>> | null {
    if (!geom)
      return;
    const extent = geom.getExtent();
    const geomFeatures = this.olMarkersSource.getFeaturesInExtent(extent)
      .filter((ft) => geom.intersectsCoordinate((ft?.getGeometry()?.getCoordinates() as Coordinate)));
      //in case of using intersectsExtent() we grab all points which contained in a "extent square"
      // .filter((ft) => ft?.getGeometry()?.intersectsExtent(extent));

    this.olSelectedMarkers.extend(geomFeatures);
    this.updateSelection(this.olSelectedMarkers);

    // return geomFeatures; //TODO: add this if we need it to handle array of selected objects
  }

  parseColorCondition(condStr: string): (any)[] {
    const resArr: (any)[] = [];
    let numbersRes = 0; //condStr.replace(/\s/g, '').match(/(\d*\.\d*|\d*)/ig);
    let symbolsRes = 0;//condStr.replace(/\s/g, '').match(/\D*/ig);
    /*
    if (!numbersRes || !symbolsRes)
      return resArr;

    numbersRes = numbersRes.filter( (e) => ((e != '') && (e != ' ')) );
    symbolsRes = symbolsRes.filter( (e) => ((e != '') && (e != ' ')) );

    if ((symbolsRes[0] === '>') || (symbolsRes[0] === '>=') ||
     (symbolsRes[0] === '<') || (symbolsRes[0] === '<=')) {
      resArr.push(symbolsRes[0]);
      resArr.push(['get', 'fieldColor']);
      resArr.push(parseFloat(numbersRes[0]));
    } else if ((symbolsRes[0] === '=') || (symbolsRes[0] === '==')) {
      resArr.push('==');
      resArr.push(['get', 'fieldColor']);
      resArr.push(parseFloat(numbersRes[0]));
    } else {
      resArr.push('between');
      resArr.push(['get', 'fieldColor']);
      resArr.push(parseFloat(numbersRes[0]));
      resArr.push(parseFloat(numbersRes[numbersRes.length-1]));
    }
    */
    return resArr;
  }

  prepareGLStyle(sizeval?: number, colorval?: string | number,
    opaval?: number, symbolval?: string, applyfilter: boolean = true): LiteralStyle {
    //
    let colorValue: any;
    let sizeValue: any;
    const minColor = toStringColor(this.markerMinColor, this.markerOpacity); //'#0000ff';
    const maxColor = toStringColor(this.markerMaxColor, this.markerOpacity); //'#ff0000';
    let colorsArray: (any)[] = [];

    if (typeof colorval !== 'undefined') {
      if (typeof colorval === 'string')
        colorValue = colorval;
      else
        colorValue = toStringColor(colorval, this.markerOpacity);
    } else if (this.useColorField === false) {
      //if color coding field doesn't assigned - we use marker defaultColor for all markers
      colorValue = toStringColor(this.defaultColor, this.markerOpacity);
    } else {
      colorValue = ['interpolate', ['linear'], ['get', 'fieldColor']]; //for linear style and for wrong colorCodingData
      //the code below create color coding for LINEAR, CONDITIONAL, CATEGORICAL color shemes from column>>
      if (this.colorCodingData !=='') {
        const parsedColors = JSON.parse(this.colorCodingData);
        if (parsedColors) {
          if (this.colorCodingType === DG.COLOR_CODING_TYPE.LINEAR) { //Linear color coding
            if (parsedColors.length > 1) {
              const interval = Math.abs(this.maxFieldColor - this.minFieldColor) / (parsedColors.length - 1);
              for (let i = 0; i < parsedColors.length; i++) {
                colorsArray.push((this.minFieldColor + (interval * i)));
                colorsArray.push(toStringColor(parsedColors[i], this.markerOpacity));
              }
            } else
              colorsArray = [this.minFieldColor, minColor, this.maxFieldColor, maxColor];
          } else if (this.colorCodingType === DG.COLOR_CODING_TYPE.CONDITIONAL) {
            //Conditional color coding
            colorValue = ['case'];
            for (const key in parsedColors) {
              if (parsedColors.hasOwnProperty(key)) {
                const parsedCondition = this.parseColorCondition(key);
                if (!parsedCondition.length)
                  continue;
                colorsArray.push(parsedCondition);
                colorsArray.push(parsedColors[key]);
              }
            }
            colorsArray.push(toStringColor(this.markerMinColor, this.markerOpacity)); //fallback value
          } else if (this.colorCodingType === DG.COLOR_CODING_TYPE.CATEGORICAL) {
            //Cathegorical color coding
            colorValue = ['match', ['get', 'fieldColor']];
            for (const key in parsedColors) {
              if (parsedColors.hasOwnProperty(key)) {
                colorsArray.push(key);
                colorsArray.push(toStringColor(parsedColors[key], this.markerOpacity));
              }
            }
            colorsArray.push(toStringColor(this.defaultColor, this.markerOpacity)); //fallback value
          }
        } else //if parsed object is wrong - we use default pattern
          colorsArray = [this.minFieldColor, minColor, this.maxFieldColor, maxColor];
      } else //if color coding string is empty - we use default pattern
        colorsArray = [this.minFieldColor, minColor, this.maxFieldColor, maxColor];

      //prepare final color coding string
      colorValue = colorValue.concat(colorsArray);
      //TODO: remove all code with old color coding above
      //new way of color coding:
      // colorValue = ['get', 'fieldColorCode']; //receive color codes stored from calls
      // colorValue = ['match', true, true, ['get', 'fieldColor'], ['get', 'fieldColor']];
      // alert(colorValue);
    }
    //OLD simple color coding approach>>
    // colorValue = ['interpolate', ['linear'], ['get', 'fieldColor'],
    //   this.minFieldColor, minColor, this.maxFieldColor, maxColor];
    //// this.minFieldColor, this.markerMinColor, this.maxFieldColor, this.markerMaxColor];

    if (typeof sizeval !== 'undefined')
      sizeValue = sizeval;
    else if (this.useSizeField === false)
      sizeValue = this.markerDefaultSize;
    else {
      sizeValue = ['interpolate', ['linear'], ['get', 'fieldSize'],
        this.minFieldSize, this.markerMinSize, this.maxFieldSize, this.markerMaxSize];
    }

    //prepare filtering condition
    let filterValue: any = true;
    if (applyfilter)
      filterValue = ['>', ['get', 'filtered'], 0];

    const markerGLStyle = {
      filter: filterValue,
      symbol: {
        symbolType: symbolval? symbolval : 'circle',
        size: sizeValue,
        color: colorValue,
        opacity: opaval? opaval : this.markerOpacity,
      },
    };
    return markerGLStyle as LiteralStyle;
  }

  updateMarkersGLLayer(recreate: boolean = true) {
    this.markerGLStyle = this.prepareGLStyle();
    this.markerGLSelStyle = this.prepareGLStyle(undefined, this.selectedColor, 0.9, undefined, false);

    // this.olMarkersLayerGL?.updateStyleVariables
    // this.olMarkersLayerGL?.setProperties({style: this.markerGLStyle});

    if (recreate) {
      let previousLayer = this.olMarkersLayerGL;
      let src = this.olMarkersLayerGL?.getSource();
      this.olMarkersLayerGL = new WebGLPointsLayer({
        source: src ? src : this.olMarkersSource,
        style: this.markerGLStyle,
      });
      this.olMarkersLayerGL.set('layerName', 'Markers GL');
      this.addLayer(this.olMarkersLayerGL);

      if (previousLayer) {
        this.olMap.removeLayer(previousLayer);
        previousLayer.dispose();
      }
      //prepare markers selection layer
      previousLayer = this.olMarkersSelLayerGL;
      src = this.olMarkersSelLayerGL?.getSource();
      this.olMarkersSelLayerGL = new WebGLPointsLayer({
        source: src ? src : this.olMarkersSelSource,
        style: this.markerGLSelStyle,
        zIndex: 100,
      });
      this.olMarkersSelLayerGL.set('layerName', 'Markers GL Selection');
      this.addLayer(this.olMarkersSelLayerGL);
      this.olMarkersSelLayerGL.setZIndex(100);

      if (previousLayer) {
        this.olMap.removeLayer(previousLayer);
        previousLayer.dispose();
      }
    } //<<if recreate
  }

  getFeatureStyleFn(feature: Feature): OLStyle.Style {
    //TODO: add here different combinations of style reading (?)
    const color = feature.get('COLOR') || '#ffeeee';
    this.styleVectorLayer.getFill().setColor(color);
    return this.styleVectorLayer;
  }

  dragNdropInteractionFn(event: any) {
    const sourceVector = new VectorSource({
      features: event.features,
    });

    OLG.addNewVectorLayer(event.file.name, null,
      function(feature) {
        const color = feature.get('COLOR') || '#ffeeee';
        OLG.styleVectorLayer.getFill().setColor(color);
        return OLG.styleVectorLayer;
      },
      sourceVector, true);
  }

  addKMLLayerFromStream(stream: string, focusOnContent: boolean = true): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource({
      features: new KML().readFeatures(stream),
    });
    const newLayer = this.addNewVectorLayer('file.name', null, null, sourceVector, focusOnContent);
    return newLayer;
  }

  addGeoJSONLayerFromStream(stream: string, focusOnContent: boolean = true): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource({
      features: new GeoJSON().readFeatures(stream),
    });
    const newLayer = this.addNewVectorLayer('file.name', null, null, sourceVector, focusOnContent);
    return newLayer;
  }

  addTopoJSONLayerFromStream(stream: string, focusOnContent: boolean = true): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource({
      features: new TopoJSON().readFeatures(stream),
    });
    const newLayer = this.addNewVectorLayer('file.name', null, null, sourceVector, focusOnContent);
    return newLayer;
  }

  addNewView(options?: Object | undefined) {
    const newView = new OLView({});
    if (options)
      newView.setProperties(options);

    this.olCurrentView = newView;
    this.olMap.setView(newView);
  }

  setViewOptions(options?: Object | undefined) {
    const oView = this.olMap.getView();
    if ((options) && (oView))
      oView.setProperties(options);
  }

  getLayersNamesList(): string[] {
    const arrayNames: string[] = [];
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        let lrName = layersArr[i].get('layerName');
        if (lrName === 'undefined')
          lrName = '';
        const layerName = lrName;
        arrayNames.push(layerName);
      }
    }
    return arrayNames;
  }

  getLayersList(): BaseLayer[] {
    return this.olMap.getAllLayers();
  }

  updateLayersList(): any[] {
    const arrLayers = this.olMap.getAllLayers();
    const arrLayersObj = [];

    for (let i = 0; i < arrLayers.length; i++) {
      const layerName = arrLayers[i].get('layerName');
      const layerId = arrLayers[i].get('layerId');
      const src = arrLayers[i].getSource();
      let exp = false;
      if (src instanceof VectorSource)
        exp = true;

      const vsb = arrLayers[i].getVisible();

      const objLayer = {
        vis: vsb,
        name: layerName,
        exp: exp,
        del: true,
        layerid: layerId,
      };

      arrLayersObj.push(objLayer);
    }

    if (this.panelLayersList)
      this.panelLayersList.refreshDF(DG.DataFrame.fromObjects(arrLayersObj));

    return arrLayersObj;
  }

  getLayerByName(lrName: string): VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|null {
    let layerResult = null;
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        if (lrName === layersArr[i].get('layerName'))
          layerResult = layersArr[i];
      }
    }
    if (layerResult instanceof WebGLPointsLayer)
      return (layerResult as WebGLPts);

    return (layerResult as VectorLayer<VectorSource>);
  }

  getLayerById(layerId: string): VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|null {
    let layerResult = null;
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        const lId = layersArr[i].get('layerId');
        if (layerId === lId)
          layerResult = layersArr[i];
      }
    }
    if (layerResult instanceof WebGLPointsLayer)
      return (layerResult as WebGLPts);

    return (layerResult as VectorLayer<VectorSource>);
  }

  removeLayerById(layerId: string): void {
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        const lId = layersArr[i].get('layerId');
        if (layerId === lId) {
          this.olMap.removeLayer(layersArr[i]);
          if (layersArr[i] instanceof WebGLPointsLayer)
            layersArr[i].dispose();
        }
      }
    }
  }

  removeAllLayers(): void {
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        this.olMap.removeLayer(layersArr[i]);
        if (layersArr[i] instanceof WebGLPointsLayer)
          layersArr[i].dispose();
      }
    }
  }

  addLayer(layerToAdd: BaseLayer) {
    layerToAdd.set('layerId', Date.now()+'|'+(Math.random()*100)); //create "UNIQUE ID"
    this.olMap.addLayer(layerToAdd);
    this.olCurrentLayer = layerToAdd;
    if (this.onRefreshCallback)
      this.onRefreshCallback();

    this.updateLayersList();
  }

  // addNewTileLayer(layerToAdd: BaseLayer) //TODO: add (if we will need it)

  //adds arbitrary Vector layer
  addNewVectorLayer(lrName?: string, opt?: Object|null,
    style?: StyleLike|null, src?: VectorSource|null, focusOnContent?: boolean): VectorLayer<VectorSource> {
    let sourceVector: VectorSource;
    if ((typeof src != 'undefined') && (src))
      sourceVector = src;
    else sourceVector = new VectorSource();
    const newLayer = new VectorLayer({source: sourceVector});

    if (lrName)
      newLayer.set('layerName', lrName);
    if (opt)
      newLayer.setProperties(opt);
    if (style)
      newLayer.setStyle(style);
    newLayer.setOpacity(this.markerOpacity);

    this.addLayer(newLayer);
    if (focusOnContent)
      this.olMap.getView().fit(sourceVector.getExtent());

    return newLayer;
  }

  //adds Bing Sattelite Map layer
  addNewBingLayer(layerName?: string | undefined, options?: Object | undefined): BaseLayer {
    const newLayer = new TileLayer({
      visible: true,
      preload: Infinity,
      source: new BingMaps({
        key: 'AkhgWhv3YTxFliztqZzt6mWy-agrbRV8EafjHeMJlCRhkIh9mwCH6k7U3hXM5e83', //TODO: hide credentials
        imagerySet: 'Aerial',
      }),
    });

    if (layerName)
      newLayer.set('layerName', layerName);
    else newLayer.set('layerName', 'Sattelite');
    if (options)
      newLayer.setProperties(options);
    this.addLayer(newLayer);
    return newLayer;
  }
  //adds Open Street Maps layer - this is our base layer
  addNewOSMLayer(layerName?: string | undefined, options?: Object | undefined): BaseLayer {
    const newLayer = new TileLayer({
      visible: true,
      preload: Infinity,
      source: new OSM()});

    if (layerName)
      newLayer.set('layerName', layerName);
    else
      newLayer.set('layerName', 'OpenStreet');
    if (options)
      newLayer.setProperties(options);
    this.addLayer(newLayer);

    return newLayer;
  }

  addNewHeatMap(layerName?: string | undefined, options?: Object | undefined): HeatmapLayer {
    const newLayer = new HeatmapLayer({
      source: this.olMarkersSource, //new VectorSource({}), //<-for now we use the same source for markers and heatmap
      blur: this.heatmapBlur,
      radius: this.heatmapRadius,
      weight: function(feature: Feature): number {
        let val = feature.get('fieldSize');
        if (typeof(val) !== 'number')
          val = 1;
        return val;
      },
    });
    if (layerName)
      newLayer.set('layerName', layerName);
    if (options)
      newLayer.setProperties(options);

    this.addLayer(newLayer);
    return newLayer;
  }

  getFeaturesFromLayer(layer: VectorLayer<VectorSource>): any[] {
    const arrFeatures: Feature[] = [];
    const src = layer.getSource();
    if (!src)
      return [];
    src.forEachFeature((ft)=>{
      arrFeatures.push(ft);
    });
    return arrFeatures;
  }

  //this function prepare text prop's - for further label placing on the map
  //we don't use it at the moment
  createTextStyle(txt: any, resolution?: number): OLStyle.Text {
    console.log('createTextStyle: ' + txt);
    return new OLStyle.Text({
      // textAlign: align == 'center',
      // font: '12px Calibri,sans-serif',
      text: txt ? (txt as string) : '',
      fill: new OLStyle.Fill({color: '#aa3300'}),
      stroke: new OLStyle.Stroke({color: '#aa3300', width: 1}),
      offsetX: 0,
      offsetY: 0,
      placement: 0,
      maxAngle: 0,
      overflow: false,
      rotation: 0,
    });
  }

  //map marker style function>>
  genStyleMarker(feature: FeatureLike, resolution: number): OLStyle.Style {
    let val = feature.get('fieldLabel');
    const size = feature.get('fieldSize');
    const clr = feature.get('fieldColor');
    if (typeof val !== 'number')
      val = 1;

    let stylePt: OLStyle.Style;
    //trying to catch renrering time and change strategy on the flight
    if (resolution < 30000 || renderTime < 1000) {
      stylePt = new OLStyle.Style({
        image: new OLStyle.Circle({
          radius: size ? size : OLG.markerDefaultSize,
          fill: new OLStyle.Fill({
            color: clr ? clr : `rgba(255, 0, 255, ${this.markerOpacity})`,
          }),
          stroke: new OLStyle.Stroke({
            color: 'rgba(255, 0, 0, 1)',
            width: 1,
          }),
        }),
        // text: 'testtext',
        // text: this.createTextStyle(val, resolution),
      });
    } else {
      stylePt = new OLStyle.Style({
        image: new OLStyle.Circle({
          radius: 1,
          stroke: new OLStyle.Stroke({
            color: 'rgba(0, 0, 0, 1)',
            width: 1,
          }),
        }),
      });
    }
    return stylePt;
  }

  static gisObjFromGeometry(ft: FeatureLike, mapref: OpenLayers | null = null):
  GisTypes.GisPoint|GisTypes.GisArea|null {
    if (!ft)
      return null;
    const geom = ft.getGeometry();
    if (geom) {
      switch (geom.getType()) {
      case 'Point': {
        const gobj = geom as OLGeom.Point;
        const coords = gobj.getCoordinates();
        const xyz = ((gobj.getLayout() == 'XYZ') || (gobj.getLayout() == 'xyz'));
        const objProps = ft.getProperties();
        if (objProps.hasOwnProperty('geometry'))
          delete objProps.geometry;
        const gisPoint = new GisTypes.GisPoint(coords[0], coords[1], xyz ? coords[2] : 0, objProps);
        return gisPoint;
      }
      case 'Polygon': {
        const gobj = geom as OLGeom.Polygon;
        const coords = gobj.getCoordinates();
        const area: GisTypes.gisPolygons = [];
        const areaPolygon: GisTypes.gisPolygon = [];
        const xyz = (gobj.getLayout().toLowerCase === 'xyz');
        //store polygon coordinates into special array>>
        for (let p = 0; p < coords.length; p++) {
          const polygonCoords: GisTypes.gisPolygonCoords = [];
          for (let i = 0; i < coords[p].length; i++) {
            const vertex = coords[p][i];
            const vertCrd: GisTypes.gisCoordinate = [vertex[0] as number, vertex[1], xyz ? vertex[2] : 0];
            polygonCoords.push(vertCrd);
          }
          areaPolygon.push(polygonCoords);
        }
        area.push(areaPolygon);
        const objProps = ft.getProperties();
        if (objProps.hasOwnProperty('geometry'))
          delete objProps.geometry;
        const gisArea = new GisTypes.GisArea(area, objProps, mapref);
        return gisArea;
      }
      case 'MultiPolygon': {
        const gobj = geom as OLGeom.MultiPolygon;
        const coords = gobj.getCoordinates();
        const area: GisTypes.gisPolygons = [];
        const xyz = (gobj.getLayout().toLowerCase === 'xyz');
        //store multipolygon coordinates into special array>>
        for (let a = 0; a < coords.length; a++) {
          const areaPolygon: GisTypes.gisPolygon = [];
          for (let p = 0; p < coords[a].length; p++) {
            const polygonCoords: GisTypes.gisPolygonCoords = [];
            for (let i = 0; i < coords[a][p].length; i++) {
              const vertex = coords[a][p][i];
              const vertCrd: GisTypes.gisCoordinate = [vertex[0] as number, vertex[1], xyz ? vertex[2] : 0];
              polygonCoords.push(vertCrd);
            }
            areaPolygon.push(polygonCoords);
          }
          area.push(areaPolygon);
        }
        const objProps = ft.getProperties();
        if (objProps.hasOwnProperty('geometry'))
          delete objProps.geometry;
        const gisArea = new GisTypes.GisArea(area, objProps, mapref);
        return gisArea;
      } //<<multi polygon converter
      } //case
    } //if geometry is valid
    return null;
  }

  exportLayerToArray(layer: VectorLayer<any>): any[] {
    const arrPreparedToDF: any[] = [];
    if (!layer)
      return arrPreparedToDF;

    const arrFeatures = this.getFeaturesFromLayer(layer);
    if (arrFeatures) {
      for (let i = 0; i < arrFeatures.length; i++) {
        const newObj = arrFeatures[i].getProperties();
        if (newObj.hasOwnProperty('geometry'))
          delete newObj.geometry;
        if (arrFeatures[i].getId())
          newObj.id_ = arrFeatures[i].getId();
        newObj.gisObject = OpenLayers.gisObjFromGeometry(arrFeatures[i], this);

        arrPreparedToDF.push(newObj);
      }
    }
    return arrPreparedToDF;
  }

  //map base events handlers>>
  onMapClick(evt: MapBrowserEvent<any>): void {
    const arrFeatures: FeatureLike[] = [];
    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
      features: arrFeatures,
    };

    this.olMap.forEachFeatureAtPixel(evt.pixel, function(feature) {
      arrFeatures.push(feature);
    });
    if (arrFeatures.length == 0)
      return;

    const ft = arrFeatures[0];
    const gisObj = OpenLayers.gisObjFromGeometry(ft, this);
    const ftGeom = ft.getGeometry();
    this.currentAreaObject = null;
    if ((ftGeom) && (ftGeom?.getType() != 'Point'))
      this.currentAreaObject = ft as Feature;
    if (gisObj) {
      // here we can invoke properties panel for our selected object (comment: it if don't need)
      setTimeout(() => {
        // grok.shell.o = DG.SemanticValue.fromValueType(gisObj, gisObj.semtype);
        grok.shell.o = gisObj;
        grok.shell.windows.showProperties = true;
      }, 50);
    }

    //USEFUL: evt.stopPropagation(); //stopImmediatePropagation();
    if (this.onClickCallback)
      this.onClickCallback(res);
  }

  onMapPointermove(evt: MapBrowserEvent<any>) {
    if (evt.dragging)
      return;

    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
      features: [],
    };

    // const ft = this.olMarkersSource.getClosestFeatureToCoordinate(evt.coordinate);
    // const ft = this.olMarkersSource.getFeaturesAtCoordinate(evt.coordinate);
    const getFeaturesOption = {
      layerFilter: this.selectCondition.bind(this),
    };
    const ftArr = this.olMap.getFeaturesAtPixel(evt.pixel, getFeaturesOption);
    if (ftArr.length > 0) {
      for (let i = 0; i < ftArr.length; i++) {
        const geom = ftArr[i].getGeometry();
        if (geom?.getType() === 'Point')
          res.features.push(ftArr[i]);
      }
    }

    if (this.onPointermoveCallback)
      this.onPointermoveCallback(res);
  }

  //map events management functions>>
  setMapSelectionCallback(f: Function) {
    this.onSelectCallback = f;
  }
  setMapClickCallback(f: Function) {
    this.onClickCallback = f;
  }
  setMapPointermoveCallback(f: Function) {
    this.onPointermoveCallback = f;
  }
  setMapRefreshCallback(f: Function) {
    this.onRefreshCallback = f;
  }
  setBtnLayersClickCallback(f: Function) {
    this.btnLayersControl.parentOnClickHandler = f;
  }

  onLayersListBtnClick() {
    if (this.panelLayersList)
      this.panelLayersList.setVisibility(this.btnLayersControl.isOn);
  }

  //map elements management functions>>
  clearLayer(layer?: VectorLayer<VectorSource> | WebGLPts | HeatmapLayer | undefined | null) {
    let aLayer: VectorLayer<VectorSource> | WebGLPts | HeatmapLayer | undefined | null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer !== 'undefined') && (layer))
      aLayer = layer;
    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.clear();
    }
  }

  addPoint(coord: Coordinate, sizeVal: number, colorVal: number,
    labelVal?: string|number|undefined, indexVal?: number|undefined,
    layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    //add marker with corresponding parameters
    let aLayer: VectorLayer<VectorSource>|HeatmapLayer|WebGLPts|undefined|null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer != 'undefined') && (layer))
      aLayer = layer;

    // const strCol = toStringColor(colorVal, this.markerOpacity);
    if (aLayer) {
      const marker = new Feature(new OLGeom.Point(OLProj.fromLonLat(coord)));
      marker.set('fieldSize', sizeVal);
      marker.set('fieldColor', colorVal);
      marker.set('fieldLabel', labelVal);
      marker.set('fieldIndex', indexVal);
      marker.set('filtered', 1);

      const src = aLayer.getSource();
      if (src)
        src.addFeature(marker);
    }
  }

  addPointFeature(feature: Feature, layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    //add marker as a predefined feature object
    if (!feature)
      return;
    let aLayer: VectorLayer<VectorSource>|HeatmapLayer|WebGLPts|undefined|null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer != 'undefined') && (layer))
      aLayer = layer;

    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.addFeature(feature);
    }
  }

  //add array of features to the layer
  addFeaturesBulk(arrFeatures: Feature<OLGeom.Geometry>[],
    layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    let aLayer: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined|null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer != 'undefined') && (layer))
      aLayer = layer;
    // const startTime = Date.now();
    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.addFeatures(arrFeatures);
    }
    // console.log('GIS addFeaturesBulk: ' + (Date.now() - startTime));
  }
}
