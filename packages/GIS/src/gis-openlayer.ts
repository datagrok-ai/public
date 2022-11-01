//base import
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as GisTypes from '../src/gis-semtypes';
import { GisViewer } from './gis-viewer';

import { Map as OLMap, MapBrowserEvent, View as OLView } from 'ol';
import HeatmapLayer from 'ol/layer/Heatmap';
import BaseLayer from 'ol/layer/Base';
import Layer from 'ol/layer/Layer';
import TileLayer from 'ol/layer/Tile';
import VectorLayer from 'ol/layer/Vector';
import WebGLPointsLayer from 'ol/layer/WebGLPoints';
import TileImage from 'ol/source/TileImage'; //this is the base class for XYZ, BingMaps etc..
import VectorSource from 'ol/source/Vector';
import Collection from 'ol/Collection';
//Projections working itilities
import * as OLProj from 'ol/proj';
import { Coordinate } from 'ol/coordinate';
//geometry drawing funtions
import Feature, {FeatureLike} from 'ol/Feature';
import * as OLGeom from 'ol/geom';
import { Type as OLType } from 'ol/geom/Geometry';
import * as OLStyle from 'ol/style';
import {LiteralStyle} from 'ol/style/literal';
//Sources import
import OSM from 'ol/source/OSM';
import BingMaps from 'ol/source/BingMaps';
import { StyleLike } from 'ol/style/Style';
//import interactions and events
import * as OLInteractions from 'ol/interaction';
import * as OLEvents from 'ol/events';
import * as OLEventsCondition from 'ol/events/condition';
import { stopPropagation } from 'ol/events/Event';

import LayerRenderer from 'ol/renderer/Layer';

//import processors
import { GPX, GeoJSON, IGC, KML, TopoJSON } from 'ol/format';
import Source from 'ol/source/Source';
import { Attribution, defaults as defaultControls } from 'ol/control';

//ZIP utilities
import JSZip, { forEach } from 'jszip';
import { zoomByDelta } from 'ol/interaction/Interaction';

export { Coordinate } from 'ol/coordinate';

type WebGLPts = WebGLPointsLayer<VectorSource<OLGeom.Point>>;

//interface for callback functions parameter
export interface OLCallbackParam {
  coord: Coordinate; //[number, number];
  pixel: [number, number];
  features: FeatureLike[];
}

let OLG: OpenLayers; //TODO: remove this terrible stuff!
let renderTime: number = 0; //temporary code for benchmarking


function toStringColor(num : number, opacity?: number) : string {
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
async function getKMLData(buffer: any): Promise<string> {
  const zip = new JSZip();
  let kmlData: string = '';
  await zip.loadAsync(buffer);
  const kmlFile = zip.file(/.kml$/i)[0];
  if (kmlFile)
    kmlData = await kmlFile.async('string');
  return kmlData;
}

export class KMZ extends KML {
  // eslint-disable-next-line camelcase
  constructor(opt_options: any) {
    // eslint-disable-next-line camelcase
    const options = opt_options || {};
    // options.iconUrlFunction = getKMLImage;
    super(options);
  }
  getType(): any {
    return 'arraybuffer'; // @typedef {'arraybuffer' | 'json' | 'text' | 'xml'} Type
  }
  readFeature(source: any, options: any) {
    const kmlData = getKMLData(source);
    return super.readFeature(kmlData, options);
  }
  readFeatures(source: any, options: any) {
    const kmlData = getKMLData(source);
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
  olMarkersSource: VectorSource<OLGeom.Point>;
  olMarkersSelSource: VectorSource<OLGeom.Point>;
  olSelectedMarkers: Collection<Feature>; //Feature<OLGeom.Point>[];

  useWebGLFlag: boolean = true;
  preventFocusing: boolean = false;
  //map interacions
  dragAndDropInteraction: OLInteractions.DragAndDrop;
  selectInteraction: OLInteractions.Select;
  dragBox: OLInteractions.DragBox;

  //styles>>
  styleVectorLayer: OLStyle.Style;
  styleVectorSelLayer: OLStyle.Style;
  //map event handlers
  onSelectCallback: Function | null = null;
  onClickCallback: Function | null = null;
  onPointermoveCallback: Function | null = null;
  onRefreshCallback: Function | null = null;

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
    this.markerGLStyle = {};
    this.markerGLSelStyle = {};
    this.olSelectedMarkers = new Collection<Feature>;

    this.olMarkersSource = new VectorSource<OLGeom.Point>();
    this.olMarkersSelSource = new VectorSource<OLGeom.Point>();

    this.styleVectorLayer = new OLStyle.Style({
      fill: new OLStyle.Fill({
        color: '#eeeeee',
        // color: 'rgba(155, 155, 55, 0.5)',
      }),
      stroke: new OLStyle.Stroke({
        // color: 'rgba(250, 250, 0, 1)',
        color: 'rgba(200, 0, 0, 1)',
        width: 1,
      }),
    });
    this.styleVectorSelLayer = new OLStyle.Style({
      fill: new OLStyle.Fill({
        color: 'rgba(200, 50, 50, 0.5)',
        // color: '#eeeeee',
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
      const extent = this.dragBox.getGeometry().getExtent();
      const boxFeatures = this.olMarkersSource.getFeaturesInExtent(extent)
        .filter((ft) => ft?.getGeometry()?.intersectsExtent(extent));

      this.olSelectedMarkers.extend(boxFeatures);
      this.updateSelection(this.olSelectedMarkers);
    });

    OLG = this;
    //end of constructor
  }

  set useWebGL(v: boolean) {
    this.useWebGLFlag = v;
    this.olMarkersLayer?.setVisible(!this.useWebGLFlag);
    this.olMarkersLayerGL?.setVisible(this.useWebGLFlag);
    this.olMarkersSelLayerGL?.setVisible(this.useWebGLFlag);
  }
  get useWebGL(): boolean { return this.useWebGLFlag; }

  set heatmapBlur(val: number) {
    this.heatmapBlurParam = val;
    if (this.olCurrentLayer instanceof HeatmapLayer)
      this.olCurrentLayer.setBlur(this.heatmapBlurParam);
    else {
      //TODO: search and apply parameter to the first heatmap layer (if this is needed)
      // this.olMap.getAllLayers()
    }
  }
  get heatmapBlur(): number { return this.heatmapBlurParam; }

  set heatmapRadius(val: number) {
    this.heatmapRadiusParam = val;
    if (this.olCurrentLayer instanceof HeatmapLayer)
      this.olCurrentLayer.setRadius(this.heatmapRadiusParam);
  }
  get heatmapRadius(): number { return this.heatmapRadiusParam; }

  selectCondition(lr: Layer<Source, LayerRenderer<any>>) {
    return ((lr === this.olMarkersLayer) || (lr === this.olMarkersLayerGL));
  }

  initMap(targetName: string) {
    if (targetName === '')
      return;

    this.olMap = new OLMap({
      target: targetName,
      controls: defaultControls({attribution: false, rotate: false}),
      view: new OLView({
        projection: 'EPSG:3857', // projection: 'EPSG:4326',
        center: OLProj.fromLonLat([-77.01072, 38.91333]),
        zoom: 5,
        enableRotation: false,
      }),
    });

    //add layers>>
    // this.addNewBingLayer('Bing sat');  //TODO: if we want to add sattelite leyer we have to manage it properly
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

    //add base event handlers>>
    this.olMap.on('click', this.onMapClick.bind(this));
    this.olMap.on('pointermove', this.onMapPointermove.bind(this));

    //add interactions to the map
    this.olMap.addInteraction(this.selectInteraction); //select interaction
    this.olMap.addInteraction(this.dragBox); //add box selection interaction
    this.olMap.addInteraction(this.dragAndDropInteraction); //add dragNdrop ability

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

  parseColorCondition(condStr: string): (any)[] {
    const resArr: (any)[] = [];
    // const regExp = /(\d*\.\d*|\d*|\D*)/ig;
    let numbersRes = condStr.replace(/\s/g, '').match(/(\d*\.\d*|\d*)/ig);
    let symbolsRes = condStr.replace(/\s/g, '').match(/\D*/ig);
    if (!numbersRes || !symbolsRes)
      return resArr;

    numbersRes = numbersRes.filter( (e) => ((e != '') && (e != ' ')) );
    symbolsRes = symbolsRes.filter( (e) => ((e != '') && (e != ' ')) );

    if ((symbolsRes[0] === '>') || (symbolsRes[0] === '>=') ||
     (symbolsRes[0] === '<') || (symbolsRes[0] === '<=')) {
      resArr.push(symbolsRes[0]);
      resArr.push(['get', 'fieldColor']);
      resArr.push(numbersRes[0]);
    } else if ((symbolsRes[0] === '=') || (symbolsRes[0] === '==')) {
      resArr.push('==');
      resArr.push(['get', 'fieldColor']);
      resArr.push(numbersRes[0]);
    } else {
      resArr.push('between');
      resArr.push(['get', 'fieldColor']);
      resArr.push(numbersRes[0]);
      resArr.push(numbersRes[numbersRes.length-1]);
    }

    return resArr;
  }

  prepareGLStyle(sizeval?: number, colorval?: string | number, opaval?: number, symbolval?: string): LiteralStyle {
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
      if (this.colorCodingData.length > 0) {
        const parsedColors = JSON.parse(this.colorCodingData);
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
      } else //if color coding string is empty - we use default pattern
        colorsArray = [this.minFieldColor, minColor, this.maxFieldColor, maxColor];

      //prepare final color coding string
      colorValue = colorValue.concat(colorsArray);
      // colorValue = ['interpolate', ['linear'], ['get', 'fieldColor'], ...colorsArray];
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

    const arrayOfFilteredRows = [];
    for (let i = 0; i < 150; i++) {
      arrayOfFilteredRows.push(i.toString());
      arrayOfFilteredRows.push(true);
    }

    const filterValue = ['match', ['get', 'fieldIndex'], ...arrayOfFilteredRows, false];

    const markerGLStyle = {
      // variables: {
      //   minVal: 0,
      // },
      // filter: filterValue,
      symbol: {
        symbolType: symbolval? symbolval : 'circle',
        size: sizeValue,
        color: colorValue,
        opacity: opaval? opaval : this.markerOpacity,
      },
    };
    return markerGLStyle as LiteralStyle;
  }

  updateMarkersGLLayer() {
    this.markerGLStyle = this.prepareGLStyle();
    this.markerGLSelStyle = this.prepareGLStyle(undefined, this.selectedColor, 0.9);

    let previousLayer = this.olMarkersLayerGL;
    let src = this.olMarkersLayerGL?.getSource();
    this.olMarkersLayerGL = new WebGLPointsLayer({
      source: src ? src : this.olMarkersSource, //new VectorSource<OLGeom.Point>(),
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
    // const tstft = src?.getFeatures();
    // alert(tstft?.length);
    this.olMarkersSelLayerGL = new WebGLPointsLayer({
      source: src ? src : this.olMarkersSelSource, //new VectorSource<OLGeom.Point>(),
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
        OLG.styleVectorLayer.getFill().setColor(color); // rgba()
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
  //adds Open Street Maps layer
  addNewOSMLayer(layerName?: string | undefined, options?: Object | undefined): BaseLayer {
    const newLayer = new TileLayer({
      visible: true,
      preload: Infinity,
      source: new OSM()});

    if (layerName)
      newLayer.set('layerName', layerName);
    else newLayer.set('layerName', 'OpenStreet');
    if (options)
      newLayer.setProperties(options);
    this.addLayer(newLayer);
    return newLayer;
  }

  addNewHeatMap(layerName?: string | undefined, options?: Object | undefined): HeatmapLayer {
    const newLayer = new HeatmapLayer({
      source: new VectorSource({}),
      blur: this.heatmapBlur,
      radius: this.heatmapRadius,
      weight: function(feature: Feature): number {
        let val = feature.get('fieldSize');
        // if (val === undefined) val = feature.get('fieldLabel');
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

  // createTextStyle(feature: FeatureLike, resolution: number): Style {
  createTextStyle(txt: any, resolution?: number): OLStyle.Text {
    console.log('createTextStyle: ' + txt);
    return new OLStyle.Text({
      // textAlign: align == 'center',
      // font: '12px Calibri,sans-serif',
      text: txt ? (txt as string) : '',
      fill: new OLStyle.Fill({ color: '#aa3300' }),
      stroke: new OLStyle.Stroke({ color: '#aa3300', width: 1 }),
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
            color: clr ? clr : `rgba(255, 0, 255, ${this.markerOpacity})`, // ${this.markerOpacity})`,
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

  static gisObjFromGeometry(ft: FeatureLike): GisTypes.GisPoint | GisTypes.GisArea | null {
    if (!ft)
      return null;
    const geom = ft.getGeometry();
    if (geom) {
      switch (geom.getType()) {
      case 'Point': {
        const gobj = geom as OLGeom.Point;
        const coords = gobj.getCoordinates();
        const xyz = ((gobj.getLayout() == 'XYZ') || (gobj.getLayout() == 'xyz'));
        const gisPoint = new GisTypes.GisPoint(coords[0], coords[1], xyz ? coords[2] : 0, ft.getProperties());
        return gisPoint;
      }
      case 'MultiPolygon':
      case 'Polygon': {
        const gobj = geom as OLGeom.Polygon;
        const coords = gobj.getCoordinates();
        const areaCoords: Array<GisTypes.gisCoordinate> = [];
        const xyz = ((gobj.getLayout() == 'XYZ') || (gobj.getLayout() == 'xyz'));
        for (let i = 0; i < coords[0].length; i++) {
          const vertex = coords[0][i]; //TODO: add inner polygons (holes) - now we use just a contour
          areaCoords.push([vertex[0] as number, vertex[1], xyz ? vertex[2] : 0]);
        }
        const gisArea = new GisTypes.GisArea(areaCoords, ft.getProperties());
        return gisArea;
      }
      //TODO: add multi polygon
      }
    }
    return null;
  }

  //map base events handlers>>
  onMapClick(evt: MapBrowserEvent<any>) {
    const arrFeatures: FeatureLike[] = [];
    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
      features: arrFeatures,
    };

    this.olMap.forEachFeatureAtPixel(evt.pixel, function(feature) {
      arrFeatures.push(feature);
    });
    if (arrFeatures.length > 0) {
      const ft = arrFeatures[0];
      const gisObj = OpenLayers.gisObjFromGeometry(ft);
      if (gisObj) { 
        //here we can invoke properties panel for our selected object (but comment it for now)
        // setTimeout(() => {
        //   // grok.shell.o = DG.SemanticValue.fromValueType(gisObj, gisObj.semtype);
        //   grok.shell.o = gisObj;
        //   grok.shell.windows.showProperties = true;
        // }, 500);
      }
    }
    // evt.stopPropagation(); //stopImmediatePropagation();
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

  //map elements management functions>>
  clearLayer(layer?: VectorLayer<VectorSource> | WebGLPts | HeatmapLayer | undefined | null) {
    let aLayer: VectorLayer<VectorSource> | WebGLPts | HeatmapLayer | undefined | null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer != 'undefined') && (layer))
      aLayer = layer;
    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.clear();
    }
  }

  addPointSc(coord: Coordinate, sizeVal: number, colorVal: number,
    labelVal?: string|number|undefined, indexVal?: number|undefined,
    layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    //add marker with corresponding parameters
    let aLayer: VectorLayer<VectorSource>|HeatmapLayer|WebGLPts|undefined|null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer != 'undefined') && (layer))
      aLayer = layer;

    const strCol = toStringColor(colorVal, this.markerOpacity);
    if (aLayer) {
      const marker = new Feature(new OLGeom.Point(OLProj.fromLonLat(coord)));
      // const style = new OLStyle.Style({
      //   image: new OLStyle.Circle({
      //     radius: sizeVal,
      //     fill: new OLStyle.Fill({
      //       color: strCol,
      //     }),
      //     stroke: new OLStyle.Stroke({
      //       color: `rgba(255, 0, 0, 1)`,
      //       width: 1,
      //     }),
      //   }),
      // });
      // marker.setStyle(style);
      marker.set('fieldSize', sizeVal);
      marker.set('fieldColor', colorVal);
      marker.set('fieldLabel', labelVal);
      marker.set('fieldIndex', indexVal);

      const src = aLayer.getSource();
      if (src)
        src.addFeature(marker);
    }
  }

  //add array of features to the layer
  addFeaturesBulk(arrFeatures: Feature<OLGeom.Geometry>[],
    layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    let aLayer: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined|null;
    aLayer = this.useWebGL ? this.olMarkersLayerGL : this.olMarkersLayer;
    if ((typeof layer != 'undefined') && (layer))
      aLayer = layer;
    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.addFeatures(arrFeatures);
    }
  }
}
