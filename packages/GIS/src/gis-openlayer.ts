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
import { useGeographic } from 'ol/proj';
import { Coordinate } from 'ol/coordinate';
//geometry drawing funtions
import * as OLGeometry from 'ol/geom/Geometry';
import { Type as OLType } from 'ol/geom/Geometry';
import * as OLPolygon from 'ol/geom/Polygon';
//@type {import("../proj/Projection.js").default|undefined}
import Feature, {FeatureLike} from 'ol/Feature';
// import * as OLFeature from 'ol/Feature';
// import Point from 'ol/geom/Point';
import * as OLGeom from 'ol/geom';
import * as OLStyle from 'ol/style';
import {LiteralStyle} from 'ol/style/literal';
//Sources import
import OSM from 'ol/source/OSM';
import BingMaps from 'ol/source/BingMaps';
import Style, { StyleLike } from 'ol/style/Style';
//import interactions and events
// import {DragAndDrop, defaults as defaultInteractions} from 'ol/interaction';
import * as OLInteractions from 'ol/interaction';
import * as OLEvents from 'ol/events';
import * as OLEventsCondition from 'ol/events/condition';

//import processors
import { GPX, GeoJSON, IGC, KML, TopoJSON } from 'ol/format';
import Source from 'ol/source/Source';
import { Attribution, defaults as defaultControls } from 'ol/control';
//ZIP utilities
import JSZip from 'jszip';
import { stopPropagation } from 'ol/events/Event';

export { Coordinate } from 'ol/coordinate';

type WebGLPts = WebGLPointsLayer<VectorSource<OLGeom.Point>>;

//interface for callback functions parameter
export interface OLCallbackParam {
  coord: Coordinate; //[number, number];
  pixel: [number, number];
}

let OLG: OpenLayers; //TODO: remove this terrible stuff!
let renderTime: number = 0;

//TODO: use this function to convert colors
//color converting
/*function hexToRGB(h: string) {
  let aRgbHex: string[] | null = ['0', '0', '0'];
  if (h.length==3) {
    aRgbHex = h.match(/.{1,1}/g);
    aRgbHex[0] = aRgbHex[0].toString().repeat(2);
    aRgbHex[1] = aRgbHex[1].toString().repeat(2);
    aRgbHex[2] = aRgbHex[2].toString().repeat(2);
  } else if (h.length==6)
    aRgbHex = h.match(/.{1,2}/g);

  const aRgb = {
    'R': parseInt(aRgbHex[0], 16),
    'G': parseInt(aRgbHex[1], 16),
    'B': parseInt(aRgbHex[2], 16),
  };
  return aRgb;
} */

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
//TODO: apply singletone pattern here?
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

  dragAndDropInteraction: OLInteractions.DragAndDrop;
  //styles>>
  styleVectorLayer: OLStyle.Style;
  styleVectorSelLayer: OLStyle.Style;
  //event handlers map
  //eventsMap = new Map<string, ()=>any>(); //TODO: solve this puzzle
  onClickCallback: Function | null = null;
  onPointermoveCallback: Function | null = null;
  onRefreshCallback: Function | null = null;
  // public labelStatus: HTMLElement | null = null;

  markerGLStyle: LiteralStyle;

  //properties from viewer
  gisViewer: GisViewer | null = null; //TODO: check do we need it for the properties access or clone prop like below
  markerDefaultSize: number = 3;
  markerMinSize: number = 1;
  markerMaxSize: number = 20;
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

  constructor(gV?: GisViewer) {
    if ((gV) && (gV instanceof GisViewer))
      this.gisViewer = gV;
    this.olMap = new OLMap({});
    this.olCurrentView = new OLView({});
    this.olCurrentLayer = null;
    this.olBaseLayer = null;
    this.olMarkersLayer = null;
    this.olMarkersLayerGL = null;
    this.markerGLStyle = {};

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
    OLG = this;
    //end of constructor
  }

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

  initMap(targetName: string) {
    if (targetName === '')
      return;

    this.olMap = new OLMap({
      target: targetName,
      controls: defaultControls({attribution: false, rotate: false}),
      view: new OLView({
        projection: 'EPSG:3857',
        // projection: 'EPSG:4326',
        center: OLProj.fromLonLat([-77.01072, 38.91333]),
        zoom: 5,
      }),
    });
    //useGeographic();

    //add layers>>
    this.addNewBingLayer('Bing sat');
    this.olBaseLayer = this.addNewOSMLayer('BaseLayer');
    // this.olMarkersLayer = this.addNewVectorLayer('Markers');
    this.olMarkersLayer = this.addNewVectorLayer('Markers', null, this.genStyleMarker.bind(this));
    this.olMarkersLayer.on('postrender', onPostrenderMarkers);
    this.olMarkersLayer.on('prerender', onPrerenderMarkers);

    this.markerGLStyle = this.prepareGLStyle();
    //TODO: this is experimental code - manage it
    this.olMarkersLayerGL = new WebGLPointsLayer({
      source: new VectorSource<OLGeom.Point>(),
      style: this.markerGLStyle,
      // opacity: this.markerOpacity,
    });
    this.olMarkersLayerGL.set('layerName', 'Markers GL');
    this.addLayer(this.olMarkersLayerGL);

    //add base event handlers>>
    this.olMap.on('click', this.onMapClick.bind(this));
    // this.olMap.on('pointermove', this.onMapPointermove);

    const selectInteraction = new OLInteractions.Select({
      condition: OLEventsCondition.click,
      style: this.styleVectorSelLayer,
    });
    this.olMap.addInteraction(selectInteraction);
    //add dragNdrop ability
    this.olMap.addInteraction(this.dragAndDropInteraction);
    //this.olMap.getInteractions().extend([selectInteraction]);
  }

  prepareGLStyle(sizeval?: number, colorval?: string, opaval?: number, symbolval?: string): LiteralStyle {
    //prepare color value:
    let colorValue: any;
    let sizeValue: any;

    if (typeof colorval !== 'undefined')
      colorValue = colorval;
    else if (this.useColorField === false)
      colorValue = toStringColor(this.defaultColor, this.markerOpacity);
    else {
      colorValue = ['interpolate', ['linear'], ['get', 'fieldColor'],
        this.minFieldColor, '#0000ff', this.maxFieldColor, '#ff0000'];
      // this.minFieldColor, this.markerMinColor, this.maxFieldColor, this.markerMaxColor];
    }

    if (typeof sizeval !== 'undefined')
      sizeValue = sizeval;
    else if (this.useSizeField === false)
      sizeValue = this.markerDefaultSize;
    else {
      sizeValue = ['interpolate', ['linear'], ['get', 'fieldSize'],
        this.minFieldSize, this.markerMinSize, this.maxFieldSize, this.markerMaxSize];
    }

    const markerGLStyle = {
      // variables: {
      //   minVal: 0,
      // },
      symbol: {
        symbolType: symbolval? symbolval : 'circle',
        // size: sizeval? sizeval : this.markerDefaultSize,
        // size: ['+', ['get', 'fieldSize'], 4],
        size: sizeValue,
        color: colorValue,
        // color: colorval? colorval : toStringColor(this.defaultColor, this.markerOpacity),
        // color: [['get', 'fieldColor']],
        opacity: opaval? opaval : this.markerOpacity,
      },
    };
    return markerGLStyle as LiteralStyle;
  }

  updateMarkersGLLayer() {
    this.markerGLStyle = this.prepareGLStyle();
    const previousLayer = this.olMarkersLayerGL;
    const src = this.olMarkersLayerGL?.getSource();
    this.olMarkersLayerGL = new WebGLPointsLayer({
      source: src ? src : new VectorSource<OLGeom.Point>(),
      style: this.markerGLStyle,
    });
    this.olMarkersLayerGL.set('layerName', 'Markers GL');
    this.olMap.addLayer(this.olMarkersLayerGL);

    if (previousLayer) {
      this.olMap.removeLayer(previousLayer);
      previousLayer.dispose();
    }
  }

  getFeatureStyleFn(feature: Feature): OLStyle.Style {
    //TODO: add here different combinations of style reading
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
      sourceVector);

    OLG.olMap.getView().fit(sourceVector.getExtent()); //TODO: check is it doubling and remove it if need
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
      //TODO: try it (for universal loader instead of different function for each of loaders)
      //var format = new ol.format.GeoJSON({
      //   defaultDataProjection: "EPSG:4326",
      //   featureProjection: "EPSG:3857"
      // });
      // var features = format.readFeatures(result);
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
          layersArr[i].dispose(); //TODO: check is it works
        }
      }
    }
  }

  removeAllLayers(): void {
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        this.olMap.removeLayer(layersArr[i]);
        // layersArr[i].unsubscribe?
        layersArr[i].dispose();
      }
    }
  }

  addLayer(layerToAdd: BaseLayer) {
    layerToAdd.set('layerId', Date.now()+'|'+(Math.random()*100));
    this.olMap.addLayer(layerToAdd);
    this.olCurrentLayer = layerToAdd;
    if (this.onRefreshCallback)
      this.onRefreshCallback();
  }

  // addNewTileLayer(layerToAdd: BaseLayer) //TODO: add

  //adds arbitrary Vector layer
  addNewVectorLayer(lrName?: string, opt?: Object|null,
    style?: StyleLike|null, src?: VectorSource|null, focusOnContent?: boolean): VectorLayer<VectorSource> {
    let sourceVector: VectorSource;
    if (src)
      sourceVector = src;
    else sourceVector = new VectorSource();
    const newLayer = new VectorLayer({source: sourceVector});

    if (lrName)
      newLayer.set('layerName', lrName);
    if (opt)
      newLayer.setProperties(opt);
    if (style)
      newLayer.setStyle(style);
    newLayer.setOpacity(0.2); //TODO: change this corresponding to layer opacity settings

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
        let val = feature.get('fieldValue');
        // if (val === undefined) val = feature.get('fieldValue');
        if (typeof(val) !== 'number')
          val = 1;
        return val;
      },
    });
    if (layerName)
      newLayer.set('layerName', layerName);
    if (options)
      newLayer.setProperties(options);
    // newLayer.changed();
    this.addLayer(newLayer);
    return newLayer;
  }

  getFeaturesFromLayer(layer: VectorLayer<VectorSource>): any[] {
    const arrFeatures: Feature[] = []; //any[]
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
  genStyleMarker(feature: FeatureLike, resolution: number): Style {
    let val = feature.get('fieldValue');
    const size = feature.get('fieldSize');
    const clr = feature.get('fieldColor');
    // if (typeof val !== 'number')
    //   val = 1;
    let stylePt: OLStyle.Style;
    //TODO: try to catch renrering time and change strategy on the flight
    if (resolution < 30000 || renderTime < 1000) {
      stylePt = new OLStyle.Style({
        image: new OLStyle.Circle({
          radius: size ? size : OLG.markerDefaultSize,
          fill: new OLStyle.Fill({
            color: clr ? clr : `rgba(255, 0, 255, 0.9)`, // ${this.markerOpacity})`,
          }),
          stroke: new OLStyle.Stroke({
            // color: 'rgba(255, 204, 0, 1)',
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
    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
    };

    const features: FeatureLike[] = [];
    OLG.olMap.forEachFeatureAtPixel(evt.pixel, function(feature) {
      features.push(feature);
    });
    if (features.length > 0) {
      const ft = features[0];
      const gisObj = OpenLayers.gisObjFromGeometry(ft);
      if (gisObj) {
        setTimeout(() => {
          // grok.shell.o = DG.SemanticValue.fromValueType(gisObj, gisObj.semtype);
          grok.shell.o = gisObj;
          grok.shell.windows.showProperties = true;
        }, 500);
      }
    }
    // evt.stopPropagation(); //stopImmediatePropagation();
    if (this.onClickCallback)
      this.onClickCallback(res);
    else {
      if (OLG) //TODO: remove this stuff (use bind)
        if (OLG.onClickCallback) OLG.onClickCallback(res); //remove this stuff (use bind)
    }
  }

  onMapPointermove(evt: MapBrowserEvent<any>) {
    if (evt.dragging)
      return;

    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
    };

    if (this.onPointermoveCallback)
      this.onPointermoveCallback(res);
    else {
      if (OLG) //TODO: remove this stuff (use bind)
        if (OLG.onPointermoveCallback) OLG.onPointermoveCallback(res); // remove this stuff (use bind)
    }
  }

  //map events management functions>>
  setMapEventHandler() {
    // this.olMap.on('click', function (evt) {
    //   displayFeatureInfo(evt.pixel);
    // });
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
    aLayer = this.olMarkersLayerGL;
    if (layer)
      aLayer = layer;
    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.clear();
    }
  }

  addPointSc(coord: Coordinate, sizeVal: number, colorVal: number,
    value?: string|number|undefined,
    layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    //
    let aLayer: VectorLayer<VectorSource>|HeatmapLayer|WebGLPts|undefined|null;
    aLayer = this.olMarkersLayer;
    if (layer)
      aLayer = layer;

    const strCol = toStringColor(colorVal, this.markerOpacity);
    if (aLayer) {
      const marker = new Feature(new OLGeom.Point(OLProj.fromLonLat(coord)));
      const style = new Style({
        image: new OLStyle.Circle({
          radius: sizeVal,
          fill: new OLStyle.Fill({
            // color: toStringColor(colorVal, this.markerOpacity), //TODO: change to gisView.markerOpacity
            color: strCol,
          }),
          stroke: new OLStyle.Stroke({
            color: `rgba(255, 0, 0, 1)`, //toStringColor(colorVal, this.markerOpacity-0.1),
            width: 1, //TODO: change to gisView.markerStrokeWidth
          }),
        }),
      });
      marker.setStyle(style);
      marker.set('fieldValue', value);

      const src = aLayer.getSource();
      if (src)
        src.addFeature(marker);
    }
  }

  addFeaturesBulk(arrFeatures: Array<Feature>,
    layer?: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined) {
    //add array of features to the layer
    let aLayer: VectorLayer<VectorSource>|WebGLPts|HeatmapLayer|undefined|null;
    aLayer = this.olMarkersLayerGL;
    if (layer)
      aLayer = layer;
    if (!aLayer)
      aLayer = this.olMarkersLayer;
    if (aLayer) {
      const src = aLayer.getSource();
      if (src)
        src.addFeatures(arrFeatures);
    }
  }
}
