import * as ui from 'datagrok-api/ui';
import * as THREE from 'three';
import ThreeGlobe from 'three-globe';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import {_package} from '../src/package.js';
import {scaleLinear, scaleSqrt, scaleSequential, interpolateYlOrRd} from 'd3';


export class GlobeViewer extends DG.JsViewer {

  constructor() {
    super();

    // Properties
    this.latitudeColumnName = this.string('latitudeColumnName');
    this.longitudeColumnName = this.string('longitudeColumnName');
    this.magnitudeColumnName = this.float('magnitudeColumnName');
    this.colorByColumnName = this.float('colorByColumnName');
    this.pointRadius = this.float('pointRadius', 15);
    this.pointAltitude = this.float('pointAltitude', 50);
    this.autorotation = this.bool('autorotation', true);

    this.rScale = scaleLinear([0, 100], [0, 1]);
    this.points = [];
    this.initialized = false;
  }

  init() {
    this.globe = new ThreeGlobe()
      .globeImageUrl(`${_package.webRoot}globe/earth-blue-marble.jpg`)
      .bumpImageUrl(`${_package.webRoot}globe/earth-topology.png`);

    this.width = this.root.parentElement.clientWidth;
    this.height = this.root.parentElement.clientHeight;

    this.renderer = new THREE.WebGLRenderer({alpha: true});
    this.renderer.domElement.style.backgroundImage = `url(${_package.webRoot}globe/night-sky.png)`;
    this.renderer.setSize(this.width, this.height);
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

    (function animate() {
      this.orbControls.update();
      this.renderer.render(this.scene, this.camera);
      requestAnimationFrame(animate.bind(this));
    }).bind(this)();

    this.initialized = true;
  }

  onTableAttached() {
    this.init();

    this.latitudeColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.LATITUDE).name;
    this.longitudeColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.LONGITUDE).name;
    this.magnitudeColumnName = this.dataFrame.columns.bySemType('Magnitude').name;
    // By default, beam color and size depend on the same column
    this.colorByColumnName = this.magnitudeColumnName;

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => {
      let width = this.root.parentElement.clientWidth;
      let height = this.root.parentElement.clientHeight;
      this.renderer.setSize(width, height);
      this.camera.aspect = width / height;
      this.camera.updateProjectionMatrix();
    }));

    this.render();
  }

  onPropertyChanged(property) {
    super.onPropertyChanged(property);
    if (this.initialized) {
      this.orbControls.autoRotate = this.autorotation;
      this.render();
    }
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  getCoordinates() {
    this.points.length = 0;
    let lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    let lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();

    let colorByCol = this.dataFrame.getCol(this.colorByColumnName);
    let mag = this.dataFrame.getCol(this.magnitudeColumnName);

    let magRange = [mag.min, mag.max];
    let color = scaleSequential().interpolator(interpolateYlOrRd);
    if (this.magnitudeColumnName === this.colorByColumnName) color.domain(magRange);
    else color.domain([colorByCol.min, colorByCol.max]);

    let altRange = [0.1, this.rScale(this.pointAltitude)];
    if (altRange[1] < 0.1) altRange[0] = 0;
    let size = scaleSqrt(magRange, altRange);

    mag = mag.getRawData();
    colorByCol = colorByCol.getRawData();
    for (let i of this.dataFrame.filter.getSelectedIndexes()) {
      this.points.push({
        lat: lat[i],
        lng: lon[i],
        size: size(mag[i]),
        color: color(colorByCol[i])
      });
    }
  }

  render() {

    this.getCoordinates();
    this.globe
      .pointsData(this.points)
      .pointAltitude('size')
      .pointColor('color')
      .pointRadius(this.rScale(this.pointRadius));

  }
}
