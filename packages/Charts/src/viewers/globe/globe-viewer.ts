import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import * as THREE from 'three';
import ThreeGlobe from 'three-globe';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import {scaleLinear, scaleSqrt, scaleSequential, interpolateYlOrRd, ScaleLinear} from 'd3';

import {_package} from '../../package';


@grok.decorators.viewer({
  name: 'Globe',
  description: 'Creates a globe viewer',
  icon: 'icons/globe-viewer.svg',
})
export class GlobeViewer extends DG.JsViewer {
  latitudeColumnName: string;
  longitudeColumnName: string;
  magnitudeColumnName: string;
  colorByColumnName: string;
  pointRadius: number;
  pointAltitude: number;
  autorotation: boolean;

  rScale: ScaleLinear<number, number, never>;
  points: any;
  initialized: boolean;

  globe?: ThreeGlobe;
  width?: number;
  height?: number;

  renderer?: THREE.WebGLRenderer;
  scene?: THREE.Scene;
  camera?: THREE.PerspectiveCamera;
  orbControls?: OrbitControls;

  constructor() {
    super();
    // Properties
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    this.magnitudeColumnName = this.string('magnitudeColumnName');
    this.colorByColumnName = this.string('colorByColumnName');
    this.pointRadius = this.float('pointRadius', 15);
    this.pointAltitude = this.float('pointAltitude', 50);
    this.autorotation = this.bool('autorotation', true);
    this.addRowSourceAndFormula();

    this.rScale = scaleLinear([0, 100], [0, 1]);
    this.points = [];
    this.initialized = false;
  }

  init() {
    this.globe = new ThreeGlobe()
      .globeImageUrl(`${_package.webRoot}img/globe/earth-blue-marble.jpg`)
      .bumpImageUrl(`${_package.webRoot}img/globe/earth-topology.png`);

    this.width = this.root.parentElement!.clientWidth;
    this.height = this.root.parentElement!.clientHeight;

    this.renderer = new THREE.WebGLRenderer({alpha: true});
    this.renderer.domElement.style.backgroundImage = `url(${_package.webRoot}img/globe/night-sky.png)`;
    this.renderer.setSize(this.width, this.height);
    if (this._testColumns())
      this.root.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.scene.add(this.globe);
    this.scene.add(new THREE.AmbientLight(0xbbbbbb));
    this.scene.add(new THREE.DirectionalLight(0xffffff, 0.6));

    this.camera = new THREE.PerspectiveCamera(50, this.width / this.height);
    this.camera.position.z = 500;

    this.orbControls = new OrbitControls(this.camera, this.renderer.domElement);
    this.orbControls.minDistance = 101;
    this.orbControls.rotateSpeed = 1.5;
    this.orbControls.zoomSpeed = 0.8;
    this.orbControls.autoRotate = true;
    this.orbControls.autoRotateSpeed = 2.2;

    (function animate(this: any) {
      this.orbControls.update();
      this.renderer.render(this.scene, this.camera);
      requestAnimationFrame(animate.bind(this));
    }).bind(this)();

    this.initialized = true;
  }

  onTableAttached() {
    this.init();
    this.filter = this.dataFrame.filter;
    this.latitudeColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.LATITUDE)?.name || '';
    this.longitudeColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.LONGITUDE)?.name || '';
    if (!this.latitudeColumnName) grok.shell.warning('Cannot find latitude column!');
    if (!this.longitudeColumnName) grok.shell.warning('Cannot find longitude column!');
    const magnitudeColumn = this.dataFrame.columns.bySemType('Magnitude')!;
    if (magnitudeColumn !== null) this.magnitudeColumnName = magnitudeColumn.name;
    else {
      const numColumns = this.dataFrame.columns.toList().filter((col) => ['double', 'int'].includes(col.type));
      if (numColumns.length !== 0)
        this.magnitudeColumnName = numColumns[0].name;
      grok.shell.info(`Cannot find magnitude column, use ${this.magnitudeColumnName} column instead`);
    }
    // By default, beam color and size depend on the same column
    this.colorByColumnName = this.magnitudeColumnName;
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => {
      const width = this.root.parentElement!.clientWidth;
      const height = this.root.parentElement!.clientHeight;
      this.renderer!.setSize(width, height);
      this.camera!.aspect = width / height;
      this.camera!.updateProjectionMatrix();
    }));
    this.render();
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    const newVal = property.get(this);
    if (property.name.endsWith('ColumnName')) {
      switch (property.name) {
      case 'latitudeColumnName':
        this.latitudeColumnName = this.dataFrame.col(newVal)?.semType === DG.SEMTYPE.LATITUDE ? newVal : '';
        break;
      case 'longitudeColumnName':
        this.longitudeColumnName = this.dataFrame.col(newVal)?.semType === DG.SEMTYPE.LONGITUDE ? newVal : '';
        break;
      case 'magnitudeColumnName':
        if (this.dataFrame.col(newVal)?.semType === 'Magnitude') this.magnitudeColumnName = newVal;
        break;
      }
    }
    if (this.initialized) {
      this.orbControls!.autoRotate = this.autorotation;
      this.render();
    }
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  getCoordinates() {
    if (!this.latitudeColumnName || !this.longitudeColumnName) {
      this.points = [];
      return;
    }
    this.points.length = 0;
    const lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    const lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();

    const colorByCol = this.dataFrame.getCol(this.colorByColumnName);
    const mag = this.dataFrame.getCol(this.magnitudeColumnName);

    const magRange = [mag.min, mag.max];
    const color = scaleSequential().interpolator(interpolateYlOrRd);
    if (this.magnitudeColumnName === this.colorByColumnName) color.domain(magRange);
    else color.domain([colorByCol.min, colorByCol.max]);

    const altRange = [0.1, this.rScale(this.pointAltitude)];
    if (altRange[1] < 0.1) altRange[0] = 0;
    const size = scaleSqrt(magRange, altRange);

    const magData = mag.getRawData();
    const colorByColData = colorByCol.getRawData();
    for (const i of this.filter.getSelectedIndexes()) {
      this.points.push({
        lat: lat[i],
        lng: lon[i],
        size: size(magData[i]),
        color: color(colorByColData[i]),
      });
    }
  }

  _testColumns() {
    const numColumns = this.dataFrame.columns.toList().filter((col) => ['double', 'int'].includes(col.type));
    return numColumns.length >= 1;
  }

  _showErrorMessage(msg: string) {this.root.appendChild(ui.divText(msg, 'd4-viewer-error'));}

  render() {
    if (!this._testColumns()) {
      this._showErrorMessage('The Globe viewer requires a minimum of 1 numerical column.');
      return;
    }

    this.getCoordinates();
    this.globe!
      .pointsData(this.points)
      .pointAltitude('size')
      .pointColor('color')
      .pointRadius(this.rScale(this.pointRadius));
  }
}
