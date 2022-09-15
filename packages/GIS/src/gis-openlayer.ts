//base import
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as GisTypes from '../src/gis-semtypes';
import {GisViewer} from './gis-viewer';

import {Map as OLMap, MapBrowserEvent, View as OLView} from 'ol';
import HeatmapLayer from 'ol/layer/Heatmap';
import BaseLayer from 'ol/layer/Base';
import Layer from 'ol/layer/Layer';
import TileLayer from 'ol/layer/Tile';
import VectorLayer from 'ol/layer/Vector';
import TileImage from 'ol/source/TileImage'; //this is the base class for XYZ, BingMaps etc..
import VectorSource from 'ol/source/Vector';
import Collection from 'ol/Collection';
//Projections working itilities
import * as OLProj from 'ol/proj';
import {useGeographic} from 'ol/proj';
import {Coordinate} from 'ol/coordinate';
//geometry drawing funtions
import * as OLGeometry from 'ol/geom/Geometry';
import {Type as OLType} from 'ol/geom/Geometry';
import * as OLPolygon from 'ol/geom/Polygon';
//@type {import("../proj/Projection.js").default|undefined}
import Feature, {FeatureLike} from 'ol/Feature';
// import * as OLFeature from 'ol/Feature';
import Point from 'ol/geom/Point';
import * as OLGeom from 'ol/geom';
import * as OLStyle from 'ol/style';
//Sources import
import OSM from 'ol/source/OSM';
import BingMaps from 'ol/source/BingMaps';
import Style, {StyleLike} from 'ol/style/Style';
//import interactions and events
// import {DragAndDrop, defaults as defaultInteractions} from 'ol/interaction';
import * as OLInteractions from 'ol/interaction';
import * as OLEvents from 'ol/events';
import * as OLEventsCondition from 'ol/events/condition';

//import processors
import {GPX, GeoJSON, IGC, KML, TopoJSON} from 'ol/format';
import Source from 'ol/source/Source';
import {Attribution, defaults as defaultControls} from 'ol/control';
//ZIP utilities
import JSZip from 'jszip';

export {Coordinate} from 'ol/coordinate';

//interface for callback functions parameter
export interface OLCallbackParam {
  coord: Coordinate; //[number, number];
  pixel: [number, number];
}

// const info = $('#info');
// info.tooltip({
//   animation: false,
//   trigger: 'manual',
// });

let OLG: OpenLayers; //TODO: remove this terrible stuff!

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

export class OpenLayers {
  olMap: OLMap;
  olCurrentView: OLView;
  olCurrentLayer: BaseLayer | VectorLayer<VectorSource> | HeatmapLayer | null;
  olBaseLayer: BaseLayer | null;
  olMarkersLayer: VectorLayer<VectorSource> | null;

  dragAndDropInteraction: OLInteractions.DragAndDrop;
  //styles>>
  styleVectorLayer: OLStyle.Style;
  styleVectorSelLayer: OLStyle.Style;
  //event handlers map
  //eventsMap = new Map<string, ()=>any>(); //TODO: solve this puzzle
  onClickCallback: Function | null = null;
  onPointermoveCallback: Function | null = null;
  // public labelStatus: HTMLElement | null = null;

  //properties from viewer
  gisViewer: GisViewer|null = null;
  markerSize: number = 1;
  markerOpacity: number = 0.7;
  weightedMarkers: boolean = false;
  heatmapBlur: number = 20;
  heatmapRadius: number = 10;

  constructor(gV?: GisViewer) {
    if ((gV) && (gV instanceof GisViewer)) this.gisViewer = gV;
    this.olMap = new OLMap({});
    this.olCurrentView = new OLView({});
    this.olCurrentLayer = null;
    this.olBaseLayer = null;
    this.olMarkersLayer = null;

    this.styleVectorLayer = new OLStyle.Style({
      fill: new OLStyle.Fill({
        // color: '#eeeeee',
        color: 'rgba(155, 155, 155, 0.5)',
      }),
      stroke: new OLStyle.Stroke({
        color: 'rgba(250, 250, 250, 0.5)',
        width: 1,
      }),
    });
    this.styleVectorSelLayer = new OLStyle.Style({
      fill: new OLStyle.Fill({
        color: 'rgba(200, 50, 50, 0.5)',
        // color: '#eeeeee',
      }),
      stroke: new OLStyle.Stroke({
        color: 'rgba(255, 0, 0, 0.7)',
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
    this.dragAndDropInteraction.on('addfeatures', this.dragNdropInteractionFn);
    OLG = this;
  }

  initMap(targetName: string) {
    if (targetName === '') return;

    this.olMap = new OLMap({
      target: targetName,
      controls: defaultControls({attribution: false, rotate: false}),
      view: new OLView({
        projection: 'EPSG:3857',
        // projection: 'EPSG:4326',
        center: OLProj.fromLonLat([34.109565, 45.452962]),
        zoom: 6,
      }),
    });
    //useGeographic();
    //add dragNdrop ability
    this.olMap.addInteraction(this.dragAndDropInteraction);

    //add layers>>
    this.addNewBingLayer('Bing sat');
    this.olBaseLayer = this.addNewOSMLayer('BaseLayer');
    this.olMarkersLayer = this.addNewVectorLayer('Markers');//, this.genStyleMarker);

    //add base event handlers>>
    this.olMap.on('click', this.onMapClick);
    // this.olMap.on('pointermove', this.onMapPointermove);

    const selectInteraction = new OLInteractions.Select({
      condition: OLEventsCondition.click,
      style: this.styleVectorSelLayer,
    });
    this.olMap.addInteraction(selectInteraction);
    //this.olMap.getInteractions().extend([selectInteraction]);
  }

  getFeatureStyleFn(feature: Feature): OLStyle.Style {
    //TODO: add here different combinations of style reading
    const color = feature.get('COLOR') || '#eeeeee';
    this.styleVectorLayer.getFill().setColor(color);
    return this.styleVectorLayer;
  }

  dragNdropInteractionFn(event: any) {
    const sourceVector = new VectorSource({
      features: event.features,
    });

    // OLG.addNewVectorLayer(event.file.name, null, null, sourceVector);
    OLG.addNewVectorLayer(event.file.name, null,
      function(feature) {
        const color = feature.get('COLOR') || '#eeeeee';
        OLG.styleVectorLayer.getFill().setColor(color); // rgba()
        return OLG.styleVectorLayer;
      },
      sourceVector);

    //TODO: if case above is workable we should add an layerID somewhere here>>
    // OLG.olMap.addLayer(
    //   new VectorLayer({
    //     source: sourceVector,
    //   }));
    // this.olMap.getView().fit(sourceVector.getExtent());
    OLG.olMap.getView().fit(sourceVector.getExtent());
  }

  addKMLLayerFromStream(stream: string, focusOnContent: boolean = true): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource({
      features: new KML().readFeatures(stream),
    });
    const newLayer = this.addNewVectorLayer('file.name', null, null, sourceVector);
    if (focusOnContent)
      this.olMap.getView().fit(sourceVector.getExtent());
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
    const newLayer = this.addNewVectorLayer('file.name', null, null, sourceVector);
    if (focusOnContent)
      this.olMap.getView().fit(sourceVector.getExtent());
    return newLayer;
  }

  addTopoJSONLayerFromStream(stream: string, focusOnContent: boolean = true): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource({
      features: new TopoJSON().readFeatures(stream),
    });
    const newLayer = this.addNewVectorLayer('file.name', null, null, sourceVector);
    if (focusOnContent)
      this.olMap.getView().fit(sourceVector.getExtent());
    return newLayer;
  }

  addNewView(options?: Object | undefined) {
    const newView = new OLView({});
    if (options) newView.setProperties(options);

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
        if (lrName === 'undefined') lrName = '';
        const layerName = lrName; //i + ' ' + lrName; //layersArr[i].get('layerName');
        arrayNames.push(layerName);
      }
    }
    return arrayNames;
  }

  getLayersList(): BaseLayer[] {
    return this.olMap.getAllLayers();
  }

  getLayerByName(lrName: string): VectorLayer<VectorSource>|HeatmapLayer|null {
    let layerResult = null;
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        if (lrName === layersArr[i].get('layerName'))
          layerResult = layersArr[i];
      }
    }
    return (layerResult as VectorLayer<VectorSource>);
  }

  getLayerById(layerId: string): VectorLayer<VectorSource>|HeatmapLayer|null {
    let layerResult = null;
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        let lId = layersArr[i].get('layerId');
        if (layerId === lId)
          layerResult = layersArr[i];
      }
    }
    return (layerResult as VectorLayer<VectorSource>);
  }

  removeLayerById(layerId: string): void {
    if (this.olMap) {
      const layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        const lId = layersArr[i].get('layerId');
        if (layerId === lId)
          this.olMap.removeLayer(layersArr[i]);
      }
    }
  }

  addLayer(layerToAdd: BaseLayer) {
    layerToAdd.set('layerId', Date.now()+'|'+(Math.random()*100));
    this.olMap.addLayer(layerToAdd);
    this.olCurrentLayer = layerToAdd;
  }

  // addNewTileLayer(layerToAdd: BaseLayer) //TODO: add

  //adds arbitrary Vector layer
  addNewVectorLayer(lrName?: string, opt?: Object|null,
    style?: StyleLike|null, src?: VectorSource|null): VectorLayer<VectorSource> {
    let sourceVector: VectorSource;
    if (src) sourceVector = src;
    else sourceVector = new VectorSource();
    const newLayer = new VectorLayer({source: sourceVector});

    if (lrName) newLayer.set('layerName', lrName);
    if (opt) newLayer.setProperties(opt);
    if (style) newLayer.setStyle(style);
    // this.olMap.addLayer(newLayer);
    newLayer.setOpacity(0.5); //TODO: change this corresponding to layer opacity settings
    this.addLayer(newLayer);
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

    if (layerName) newLayer.set('layerName', layerName);
    else newLayer.set('layerName', 'Sattelite');
    if (options) newLayer.setProperties(options);
    // this.olMap.addLayer(newLayer);
    this.addLayer(newLayer);
    return newLayer;
  }
  //adds Open Street Maps layer
  addNewOSMLayer(layerName?: string | undefined, options?: Object | undefined): BaseLayer {
    const newLayer = new TileLayer({
      visible: true,
      preload: Infinity,
      source: new OSM()});

    if (layerName) newLayer.set('layerName', layerName);
    else newLayer.set('layerName', 'OpenStreet');
    if (options) newLayer.setProperties(options);
    // this.olMap.addLayer(newLayer);
    this.addLayer(newLayer);
    return newLayer;
  }

  addNewHeatMap(layerName?: string | undefined, options?: Object | undefined): HeatmapLayer {
    const newLayer = new HeatmapLayer({
      source: new VectorSource({}),
      blur: this.gisViewer ? this.gisViewer.heatmapBlur : this.heatmapBlur,
      radius: this.gisViewer ? this.gisViewer.heatmapBlur : this.heatmapRadius,
      weight: function(feature: Feature): number {
        let val = feature.get('fieldValue');
        if (typeof(val) !== 'number') val = 1;
        return val;
      },
    });
    if (layerName) newLayer.set('layerName', layerName);
    if (options) newLayer.setProperties(options);
    // newLayer.changed();
    this.addLayer(newLayer);
    return newLayer;
  }

  getFeaturesFromLayer(layer: VectorLayer<VectorSource>): any[] {
    const featuresColl: any[] = [];
    const src = layer.getSource();
    if (!src) return [];
    src.forEachFeature((ft)=>{
      featuresColl.push(ft.getProperties());
    });
    return featuresColl;
  }

  //map marker style function>>
  genStyleMarker(feature: Feature): Style {
    let val = feature.get('fieldValue');
    if (typeof(val) !== 'number') val = 1;

    const style = new Style({
      image: new OLStyle.Circle({
        radius: this.weightedMarkers ? val*1 : this.markerSize,
        fill: new OLStyle.Fill({
          color: 'rgba(255, 153, 0, 0.4)',
        }),
        stroke: new OLStyle.Stroke({
          color: 'rgba(255, 204, 0, 0.2)',
          width: 1,
        }),
      }),
    });
    return style;
  }

  //map base events handlers>>
  onMapClick(evt: MapBrowserEvent<any>) {
    // evt.coordinate
    //evt.pixel
    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
    };

    const features: FeatureLike[] = [];
    OLG.olMap.forEachFeatureAtPixel(evt.pixel, function(feature) {
      features.push(feature);
    });
    if (features.length > 0) {
      let ft = features[0];
      const geom = ft.getGeometry();
      if (geom) {
        switch (geom.getType()) {
        case 'Point': {
          const gobj = geom as OLGeom.Point;
          const coords = gobj.getCoordinates();
          const xyz = ((gobj.getLayout() == 'XYZ') || (gobj.getLayout() == 'xyz'));
          const gisPoint = new GisTypes.GisPoint(coords[0], coords[1], xyz ? coords[2] : 0, ft.getProperties());
          grok.shell.o = gisPoint;
          break;
        }
        case 'Polygon': {
          const gobj = geom as OLGeom.Polygon;
          const coords = gobj.getCoordinates();
          const areaCoords: Array<GisTypes.gisCoordinate> = [];
          const xyz = ((gobj.getLayout() == 'XYZ') || (gobj.getLayout() == 'xyz'));
          for (let i = 0; i < coords[0].length; i++) {
            let vertex = coords[0][i]; //TODO: add inner polygons (holes) - now we use just a contour
            areaCoords.push([vertex[0] as number, vertex[1], xyz ? vertex[2] : 0]);
          }
          const gisArea = new GisTypes.GisArea(areaCoords, ft.getProperties());
          grok.shell.o = gisArea;
          break;
        }
        }
      }
    }
    if (this.onClickCallback)
      this.onClickCallback(res);
    else {
      if (OLG) //TODO: remove this stuff (use bind)
        if (OLG.onClickCallback) OLG.onClickCallback(res);
    }
  }

  onMapPointermove(evt: MapBrowserEvent<any>) {
    if (evt.dragging) return;

    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]],
    };

    if (this.onPointermoveCallback)
      this.onPointermoveCallback(res);
    else {
      if (OLG) //TODO: remove this stuff (use bind)
        if (OLG.onPointermoveCallback) OLG.onPointermoveCallback(res);
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

  //map elements management functions>>
  clearLayer(layer?: VectorLayer<VectorSource> | HeatmapLayer | undefined | null) {
    let aLayer: VectorLayer<VectorSource> | HeatmapLayer | undefined | null;
    aLayer = this.olMarkersLayer;
    if (layer) aLayer = layer;
    if (aLayer) {
      const src = aLayer.getSource();
      if (src) src.clear();
    }
  }

  addPoint(coord: Coordinate, value?: string|number|undefined,
    layer?: VectorLayer<VectorSource>|HeatmapLayer|undefined) {
    //
    let aLayer: VectorLayer<VectorSource>|HeatmapLayer|undefined|null;
    aLayer = this.olMarkersLayer;
    if (layer) aLayer = layer;

    const maxradius = 30;
    if (aLayer) {
      let val = value;
      if (typeof(val) !== 'number') val = 1;
      let rad = this.markerSize;
      if (this.weightedMarkers) {
        rad = val;
        if (val > maxradius) rad = maxradius;
      }
      let stroke = 1;
      if (this.weightedMarkers) {
        // (val>maxradius) ? ((val/10>maxradius/2) ? (maxradius/2+2) : (val/10) ) : 1
        stroke = val;
        if ((val > maxradius) && (val/10 > maxradius/2)) stroke = maxradius/2+2;
        else stroke = val/10;
      }
      const marker = new Feature(new Point(OLProj.fromLonLat(coord)));
      const style = new Style({
        image: new OLStyle.Circle({
          radius: rad,
          fill: new OLStyle.Fill({
            color: `rgba(255, 153, 0, ${this.markerOpacity})`,
          }),
          stroke: new OLStyle.Stroke({
            color: `rgba(255, 204, 0, ${this.markerOpacity-0.2})`,
            width: stroke,
            // width: (val>maxradius) ? ((val/10>maxradius/2) ? (maxradius/2+2) : (val/10) ) : 1,
          }),
        }),
      });
      marker.setStyle(style);
      marker.set('fieldValue', val);

      const src = aLayer.getSource();
      if (src) src.addFeature(marker);
    }
  }

  addPointSc(coord: Coordinate, sizeVal: number, colorVal: number,
    value?: string|number|undefined,
    layer?: VectorLayer<VectorSource>|HeatmapLayer|undefined) {
    //
    let aLayer: VectorLayer<VectorSource>|HeatmapLayer|undefined|null;
    aLayer = this.olMarkersLayer;
    if (layer) aLayer = layer;

    if (aLayer) {
      const marker = new Feature(new Point(OLProj.fromLonLat(coord)));
      const style = new Style({
        image: new OLStyle.Circle({
          radius: sizeVal,
          fill: new OLStyle.Fill({
            color: toStringColor(colorVal, this.markerOpacity), //TODO: change to gisView.markerOpacity
          }),
          stroke: new OLStyle.Stroke({
            color: `rgba(255, 204, 0, 0.4)`, //toStringColor(colorVal, this.markerOpacity-0.1),
            width: 1, //TODO: change to gisView.markerStrokeWidth
          }),
        }),
      });
      marker.setStyle(style);
      marker.set('fieldValue', value);

      const src = aLayer.getSource();
      if (src) src.addFeature(marker);
    }
  }
}
