import {Map, View as OLView} from 'ol';
import TileLayer from 'ol/layer/Tile';
import TileImage from 'ol/source/TileImage'; //this is the base class for XYZ, BingMaps etc..
import VectorLayer from 'ol/layer/Vector';
import VectorSource from 'ol/source/Vector';
//Projections working itilities
import * as OLProj from 'ol/proj';
import {Coordinate} from 'ol/coordinate';
//geometry drawing funtions
import * as OLPolygon from 'ol/geom/Polygon';
import Feature from 'ol/Feature';
import Point from 'ol/geom/Point';
//Sources import
import OSM from 'ol/source/OSM';
import BingMaps from 'ol/source/BingMaps';
import BaseLayer from 'ol/layer/Base';

export class OpenLayers {
  olMap: Map;
  olLayers: TileLayer<TileImage>[]; //array or map layers TODO: change it to Collection ??
  olCurrentView: OLView;
  olCurrentLayer: BaseLayer | VectorLayer<VectorSource> | null;
  olBaseLayer: BaseLayer | null;
  olMarkersLayer: VectorLayer<VectorSource> | null;
  constructor() {
    this.olMap = new Map({});
    this.olLayers = [];
    this.olCurrentView = new OLView({});
    this.olCurrentLayer = null;
    this.olBaseLayer = null;
    this.olMarkersLayer = null;
  }
  initMap(targetName: string) {
    if (targetName === '') return;

    this.olMap = new Map({
      target: targetName,
      //layers: this.olLayers,
      view: new OLView({
        center: OLProj.fromLonLat([34.109565029604006, 45.45296242157734]),
        zoom: 7,
      }),
    });

    // this.addNewView({
    //   center: OLProj.fromLonLat([34.109565029604006, 45.45296242157734]),
    //   zoom: 7,
    // });
    this.olBaseLayer = this.addNewOSMLayer('BaseLayer');
    this.olMarkersLayer = this.addNewVectorLayer('Markers');
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

  // addNewTileLayer(layerToAdd: BaseLayer) //TODO: add

  //adds arbitrary Vector layer
  addNewVectorLayer(layerName?: string | undefined, options?: Object | undefined): VectorLayer<VectorSource> {
    const sourceVector = new VectorSource();
    const newLayer = new VectorLayer({source: sourceVector});

    if (layerName) newLayer.set('layerName', layerName);
    if (options) newLayer.setProperties(options);
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

  //map elements management functions>>
  addPoint(coord: Coordinate, layer?: VectorLayer<VectorSource> | undefined) {
    let aLayer = this.olMarkersLayer;
    if (layer) aLayer = layer;

    if (aLayer) {
      const src = aLayer.getSource();
      src?.addFeature(new Feature(new Point(OLProj.fromLonLat(coord))));
    }
  }
}
