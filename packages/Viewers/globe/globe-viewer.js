import * as THREE from 'three';
import ThreeGlobe from 'three-globe';
import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls.js';


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
        for (let i = 0; i < rowCount; i++) {
            points.push({
                lat: latCol[i],
                lng: lonCol[i],
                size: magCol[i]
            });
        }
        let globe = new ThreeGlobe()
            .globeImageUrl('//unpkg.com/three-globe/example/img/earth-dark.jpg')
            .bumpImageUrl('//unpkg.com/three-globe/example/img/earth-topology.png')
            .pointsData(points)
            .pointAltitude('size');

        // Basic example
        let width = this.root.parentElement.clientWidth;
        let height = this.root.parentElement.clientHeight;
        
        // Setup renderer
        let renderer = new THREE.WebGLRenderer();
        renderer.setSize(width, height);
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
        let tbControls = new TrackballControls(camera, renderer.domElement);
        tbControls.minDistance = 101;
        tbControls.rotateSpeed = 5;
        tbControls.zoomSpeed = 0.8;

        // Kick-off renderer
        (function animate() { // IIFE
            // Frame cycle
            tbControls.update();
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
