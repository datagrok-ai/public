import {Map as OLMap, MapBrowserEvent, View as OLView} from 'ol';
import HeatmapLayer from 'ol/layer/Heatmap';
import BaseLayer from 'ol/layer/Base';
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
import KML from 'ol/format/KML';
import Style, {StyleLike} from 'ol/style/Style';

export {Coordinate} from 'ol/coordinate';
//interface for callback functions parameter
export interface OLCallbackParam {
  coord: Coordinate; //[number, number];
  pixel: [number, number];
}

export class OpenLayers {
  olMap: OLMap;
  olCurrentView: OLView;
  olCurrentLayer: BaseLayer | VectorLayer<VectorSource> | HeatmapLayer | null;
  olBaseLayer: BaseLayer | null;
  olMarkersLayer: VectorLayer<VectorSource> | null;

  //event handlers map
  eventsMap = new Map<string, ()=>any>(); //TODO: solve this puzzle
  onClickCallback: Function | null = null;
  onPointermoveCallback: Function | null = null;
  public labelStatus: HTMLElement | null = null;

  constructor() {
    this.olMap = new OLMap({});
    this.olCurrentView = new OLView({});
    this.olCurrentLayer = null;
    this.olBaseLayer = null;
    this.olMarkersLayer = null;
  }

  initMap(targetName: string) {
    if (targetName === '') return;

    this.olMap = new OLMap({
      target: targetName,
      view: new OLView({
        center: OLProj.fromLonLat([34.109565, 45.452962]),
        zoom: 7,
      }),
    });

    //add layers>>
    this.olBaseLayer = this.addNewOSMLayer('BaseLayer');
    this.olMarkersLayer = this.addNewVectorLayer('Markers');//, this.genStyleMarker);

    //add base event handlers>>
    this.olMap.on('click', this.onMapClick);
    this.olMap.on('pointermove', this.onMapPointermove);
  }

  addNewView(options?: Object | undefined) {
    const newView = new OLView({});
    if (options) newView.setProperties(options);

    this.olCurrentView = newView;
    this.olMap.setView(newView);
  }

  addLayer(layerToAdd: BaseLayer) {
    this.olMap.addLayer(layerToAdd);
    this.olCurrentLayer = layerToAdd;
  }

  getLayersList(): string[] {
    const arrayNames: string[] = [];

    if (this.olMap) {
      let layersArr = this.olMap.getAllLayers();
      for (let i = 0; i < layersArr.length; i++) {
        let layerName = i + ' ' + layersArr[i].get('layerName');
        //if (!layerName) layerName = i + ':';
        arrayNames.push(layerName);
      }
    }
    return arrayNames;
  }

  // addNewTileLayer(layerToAdd: BaseLayer) //TODO: add

  //adds arbitrary Vector layer
  addNewVectorLayer(lrName?: string, opt?: Object, style?: StyleLike): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource();
    const newLayer = new VectorLayer({source: sourceVector});

    if (lrName) newLayer.set('layerName', lrName);
    if (opt) newLayer.setProperties(opt);
    if (style) newLayer.setStyle(style);
    this.olMap.addLayer(newLayer);
    return newLayer;
  }

  //adds Open Street Maps layer
  addNewOSMLayer(layerName?: string | undefined, options?: Object | undefined): BaseLayer {
    const newLayer = new TileLayer({
      visible: true,
      preload: Infinity,
      source: new OSM()});

    if (layerName) newLayer.set('layerName', layerName);
    if (options) newLayer.setProperties(options);
    this.olMap.addLayer(newLayer);
    return newLayer;
  }

  addNewHeatMap(layerName?: string | undefined, options?: Object | undefined): HeatmapLayer {
    const newLayer = new HeatmapLayer({
      source: new VectorSource({}),
      blur: 25,
      radius: 10,
      weight: function(feature: Feature): number {
        let val = feature.get('fieldValue');
        if (typeof(val) !== 'number') val = 1;
        return val;
      },
    });
    if (layerName) newLayer.set('layerName', layerName);
    if (options) newLayer.setProperties(options);
    this.olMap.addLayer(newLayer);
    return newLayer;
  }

  //map marker style function>>
  genStyleMarker(feature: Feature): Style {
    let val = feature.get('fieldValue');
    if (typeof(val) !== 'number') val = 1;

    const style = new Style({
      image: new OLStyle.Circle({
        radius: val*20,
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
      pixel: [evt.pixel[0], evt.pixel[1]]
    };

    if (this.onClickCallback)
      this.onClickCallback(res);
  }
  onMapPointermove(evt: MapBrowserEvent<any>) {
    if (evt.dragging) return;

    const res: OLCallbackParam = {
      coord: evt.coordinate,
      pixel: [evt.pixel[0], evt.pixel[1]]
    };

    //TODO: remove this crutch - only callback fn
    if (this.labelStatus)
      this.labelStatus.innerText = evt.coordinate[0] + ', ' + evt.coordinate[1];
    else {
      //TODO: remove this crutch
      let lbl = document.getElementById('lbl-coord');
      if (lbl) lbl.innerHTML = evt.coordinate[0] + ', ' + evt.coordinate[1];
    }

    if (this.onPointermoveCallback)
      this.onPointermoveCallback(res);
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
  clearLayer(layer?: VectorLayer<VectorSource> | HeatmapLayer | undefined) {
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

    if (aLayer) {
      let val = value;
      if (typeof(val) !== 'number') val = 1;

      const marker = new Feature(new Point(OLProj.fromLonLat(coord)));
      const style = new Style({
        image: new OLStyle.Circle({
          radius: val*1,
          fill: new OLStyle.Fill({
            color: 'rgba(255, 153, 0, 0.4)',
          }),
          stroke: new OLStyle.Stroke({
            color: 'rgba(255, 204, 0, 0.2)',
            width: 1,
          }),
        }),
      });
      marker.setStyle(style);

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
