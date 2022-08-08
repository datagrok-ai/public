import {Map as OLMap, MapBrowserEvent, View as OLView} from 'ol';
import HeatmapLayer from 'ol/layer/Heatmap';
import BaseLayer from 'ol/layer/Base';
import Layer from 'ol/layer/Layer';
import TileLayer from 'ol/layer/Tile';
import VectorLayer from 'ol/layer/Vector';
import TileImage from 'ol/source/TileImage'; //this is the base class for XYZ, BingMaps etc..
import VectorSource from 'ol/source/Vector';
//Projections working itilities
import * as OLProj from 'ol/proj';
import {Coordinate} from 'ol/coordinate';
//geometry drawing funtions
import * as OLPolygon from 'ol/geom/Polygon';
import Feature from 'ol/Feature';
import Point from 'ol/geom/Point';
import * as OLStyle from 'ol/style';
//Sources import
import OSM from 'ol/source/OSM';
import BingMaps from 'ol/source/BingMaps';
import Style, {StyleLike} from 'ol/style/Style';
//import processors
import {DragAndDrop, defaults as defaultInteractions} from 'ol/interaction';
import {GPX, GeoJSON, IGC, KML, TopoJSON} from 'ol/format';
import Source from 'ol/source/Source';
import {Attribution, defaults as defaultControls} from 'ol/control';

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

let OLG: OpenLayers; //TODo: remove this terrible stuff!

export class OpenLayers {
  olMap: OLMap;
  olCurrentView: OLView;
  olCurrentLayer: BaseLayer | VectorLayer<VectorSource> | HeatmapLayer | null;
  olBaseLayer: BaseLayer | null;
  olMarkersLayer: VectorLayer<VectorSource> | null;

  dragAndDropInteraction: DragAndDrop;
  //event handlers map
  //eventsMap = new Map<string, ()=>any>(); //TODO: solve this puzzle
  onClickCallback: Function | null = null;
  onPointermoveCallback: Function | null = null;
  // public labelStatus: HTMLElement | null = null;

  //properties from viewer
  markerSize: number = 1;
  markerOpacity: number = 0.7;
  weightedMarkers: boolean = false;
  heatmapBlur: number = 20;
  heatmapRadius: number = 10;

  constructor() {
    this.olMap = new OLMap({});
    this.olCurrentView = new OLView({});
    this.olCurrentLayer = null;
    this.olBaseLayer = null;
    this.olMarkersLayer = null;

    this.dragAndDropInteraction = new DragAndDrop({
      formatConstructors: [
        KML,
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
        center: OLProj.fromLonLat([34.109565, 45.452962]),
        zoom: 7,
      }),
    });

    //add dragNdrop ability
    this.olMap.addInteraction(this.dragAndDropInteraction);

    //add layers>>
    this.addNewBingLayer('Bing sat');
    this.olBaseLayer = this.addNewOSMLayer('BaseLayer');
    this.olMarkersLayer = this.addNewVectorLayer('Markers');//, this.genStyleMarker);

    //add base event handlers>>
    this.olMap.on('click', this.onMapClick);
    this.olMap.on('pointermove', this.onMapPointermove);
  }

  dragNdropInteractionFn(event: any) {
    const sourceVector = new VectorSource({
      features: event.features,
    });
    OLG.addNewVectorLayer(event.file.name, null, null, sourceVector);
    //TODO: if case above is workable we should add an layerID somewhere here>>
    // OLG.olMap.addLayer(
    //   new VectorLayer({
    //     source: sourceVector,
    //   }));
    // this.olMap.getView().fit(sourceVector.getExtent());
    OLG.olMap.getView().fit(sourceVector.getExtent());
  }

  addKMLLayerFromStream(stream: string) {
    const sourceVector = new VectorSource({
      features: new KML().readFeatures(stream),
    });
    this.addNewVectorLayer('file.name', null, null, sourceVector);
    this.olMap.getView().fit(sourceVector.getExtent());
  }

  addNewView(options?: Object | undefined) {
    const newView = new OLView({});
    if (options) newView.setProperties(options);

    this.olCurrentView = newView;
    this.olMap.setView(newView);
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
    this.addLayer(newLayer);
    return newLayer;
  }

  //adds Bing Sattelite Map layer
  addNewBingLayer(layerName?: string | undefined, options?: Object | undefined): BaseLayer {
    const newLayer = new TileLayer({
      visible: true,
      preload: Infinity,
      source: new BingMaps({
        key: 'AkhgWhv3YTxFliztqZzt6mWy-agrbRV8EafjHeMJlCRhkIh9mwCH6k7U3hXM5e83',
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
      blur: this.heatmapBlur,
      radius: this.heatmapRadius,
      weight: function(feature: Feature): number {
        let val = feature.get('fieldValue');
        if (typeof(val) !== 'number') val = 1;
        return val;
      },
    });
    if (layerName) newLayer.set('layerName', layerName);
    if (options) newLayer.setProperties(options);
    // this.olMap.addLayer(newLayer);
    this.addLayer(newLayer);
    return newLayer;
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

  // displayFeatureInfo = function (pixel) {
  //   info.css({
  //     left: pixel[0] + 'px',
  //     top: pixel[1] - 15 + 'px',
  //   });
  //   const feature = map.forEachFeatureAtPixel(pixel, function (feature) {
  //     return feature;
  //   });
  //   if (feature) {
  //     info.attr('data-original-title', feature.get('name')).tooltip('show');
  //   } else {
  //     info.tooltip('hide');
  //   }
  // };

  //map base events handlers>>
  onMapClick(evt: MapBrowserEvent<any>) {
    // evt.coordinate
    //evt.pixel
    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]]
    };

    const features = this.olMap.getFeaturesAtPixel(evt.pixel);
    if (features) {
      //(features[0].getGeometry).
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
      pixel: [evt.pixel[0], evt.pixel[1]]
    };

    // const feature = this.olMap.forEachFeatureAtPixel(evt.pixel, function(feature) {
    //   return feature;
    // });
    // let featuredata: string = ' / Value';
    // if (feature) featuredata = featuredata + feature.get('fieldValue');
    //TODO: remove this crutch - only callback fn
    // let lbl = document.getElementById('lbl-coord');
    // if (lbl) lbl.innerHTML = evt.coordinate[0] + ', ' + evt.coordinate[1];

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

    const maxradius = 60;
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
  //   addPointsSet(coord: Coordinate, layer?: VectorLayer<VectorSource> | undefined) {
  //     let aLayer = this.olMarkersLayer;
  //     if (layer) aLayer = layer;

//     if (aLayer) {
//       const src = aLayer.getSource();
// //      src?.addFeatures()
//       src?.addFeature(new Feature(new Point(OLProj.fromLonLat(coord))));
//     }
//   }
}
