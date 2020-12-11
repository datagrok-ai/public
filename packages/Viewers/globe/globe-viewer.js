import * as THREE from 'three';
import ThreeGlobe from 'three-globe';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { _package } from '../src/package.js';


export class GlobeViewer extends DG.JsViewer {

    constructor() {
        super();

        // Properties
        this.latitude = this.string('latitudeColumnName');
        this.longitude = this.string('longitudeColumnName');
        this.magnitude = this.float('magnitudeColumnName');
    }

    init() {

        let latCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.LATITUDE).getRawData();
        let lonCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.LONGITUDE).getRawData();
        let magCol = this.dataFrame.columns.toList().filter(col => col.type === 'double')[0].getRawData();

        let points = [];
        let rowCount = this.dataFrame.rowCount;
        let factor = 0.005;
        for (let i = 0; i < rowCount; i++) {
            points.push({
                lat: latCol[i],
                lng: lonCol[i],
                size: magCol[i] * factor,
                color: '#ffd444'
            });
        }
        let globe = new ThreeGlobe()
            .globeImageUrl(`${_package.webRoot}globe/earth-blue-marble.jpg`)
            .bumpImageUrl(`${_package.webRoot}globe/earth-topology.png`)
            .pointsData(points)
            .pointAltitude('size')
            .pointColor('color')
            .pointRadius(0.15);

        // Basic example
        let width = this.root.parentElement.clientWidth;
        let height = this.root.parentElement.clientHeight;
        
        // Setup renderer
        let renderer = new THREE.WebGLRenderer({ alpha: true });
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
        camera.aspect = width/height;
        camera.updateProjectionMatrix();
        camera.position.z = 500;

        // Add camera controls
        let orbControls = new OrbitControls(camera, renderer.domElement);
        orbControls.autoRotate = true;
        orbControls.autoRotateSpeed = 2.2;
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

    onTableAttached() {
        this.init();
        this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
        this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));
        this.render();
    }

    render() {}
}
