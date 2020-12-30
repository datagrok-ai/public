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
    // TODO: draw country polygons for columns with DG.SEMTYPE.COUNTRY
    this.magnitudeColumnName = this.float('magnitudeColumnName');
    this.pointRadius = this.float('pointRadius', 15);
    this.pointAltitude = this.float('pointAltitude', 50);
    this.autorotation = this.bool('autorotation', true);

    this.rScale = scaleLinear([0, 100], [0, 1]);
    this.points = [];
    this.initialized = false;
  }

  init() {
    this.initialized = true;
  }

  onTableAttached() {
    this.init();

    this.latitudeColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.LATITUDE).name;
    this.longitudeColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.LONGITUDE).name;
    this.magnitudeColumnName = this.dataFrame.columns.bySemType('Magnitude').name;
    this.getCoordinates();

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));

    this.render();
  }

  onPropertyChanged(property) {
    super.onPropertyChanged(property);
    if (this.initialized) {
      this.getCoordinates();
      this.render();
    }
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  getCoordinates() {
    let lat = this.dataFrame.getCol(this.latitudeColumnName).getRawData();
    let lon = this.dataFrame.getCol(this.longitudeColumnName).getRawData();
    let mag = this.dataFrame.getCol(this.magnitudeColumnName);
    let magRange = [mag.min, mag.max];
    let color = scaleSequential(magRange, interpolateYlOrRd);
    let altRange = [0.1, this.rScale(this.pointAltitude)];
    if (altRange[1] < 0.1) altRange[0] = 0;
    let size = scaleSqrt(magRange, altRange);
    mag = mag.getRawData();
    let rowCount = this.dataFrame.rowCount;
    for (let i = 0; i < rowCount; i++) {
      this.points.push({
        lat: lat[i],
        lng: lon[i],
        size: size(mag[i]),
        color: color(mag[i])
      });
    }
  }

  render() {

    let globe = new ThreeGlobe()
      .globeImageUrl(`${_package.webRoot}globe/earth-blue-marble.jpg`)
      .bumpImageUrl(`${_package.webRoot}globe/earth-topology.png`)
      .pointsData(this.points)
      .pointAltitude('size')
      .pointColor('color')
      .pointRadius(this.rScale(this.pointRadius));

    $(this.root).empty();
    let width = this.root.parentElement.clientWidth;
    let height = this.root.parentElement.clientHeight;

    // Setup renderer
    let renderer = new THREE.WebGLRenderer({alpha: true});
    renderer.setSize(width, height);
    renderer.domElement.style.backgroundImage = `url(${_package.webRoot}globe/night-sky.png)`;
    this.root.appendChild(renderer.domElement);

    // Setup scene
    let scene = new THREE.Scene();
    scene.add(globe);
    scene.add(new THREE.AmbientLight(0xbbbbbb));
    scene.add(new THREE.DirectionalLight(0xffffff, 0.6));

    // Setup camera
    let camera = new THREE.PerspectiveCamera();
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    camera.position.z = 500;

    // Add camera controls
    let orbControls = new OrbitControls(camera, renderer.domElement);
    if (this.autorotation) {
      orbControls.autoRotate = true;
      orbControls.autoRotateSpeed = 2.2;
    }
    orbControls.minDistance = 101;
    orbControls.rotateSpeed = 1.5;
    orbControls.zoomSpeed = 0.8;

    // Kick-off renderer
    (function animate() { // IIFE
      // Frame cycle
      orbControls.update();
      renderer.render(scene, camera);
      requestAnimationFrame(animate);
    })();
  }
}
