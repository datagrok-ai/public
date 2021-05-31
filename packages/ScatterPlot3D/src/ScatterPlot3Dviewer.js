import * as THREE from 'three'
//import { OrbitControls } from "https://threejs.org/examples/jsm/controls/OrbitControls.js";
 
export class ScatterPlot3Dviewer extends DG.JsViewer {
    constructor() {
      super();
        console.log("HREEE ", THREE)
      // Register properties and define fields initialized to properties' default values
      // Properties that represent columns should end with the 'ColumnName' postfix
      this.splitColumnName = this.string('splitColumnName', 'site');
      this.valueColumnName = this.int('valueColumnName', 'age');
  
      // Provide a list of acceptable values in the field `choices`
      this.valueAggrType = this.string('valueAggrType', 'avg', { choices: ['avg', 'count', 'sum'] });
      this.color = this.string('color', 'steelblue', { choices: ['darkcyan', 'seagreen', 'steelblue'] });
      this.initialized = false;
	  this.initLayout();
    }
	
	rendererResize(size) {
		this.renderer.setSize(size.width, size.height);
	}

    initLayout() {
		this.camera = new THREE.PerspectiveCamera(45, 1, 2, 100000);
		this.scene = new THREE.Scene();
		this.camera.position.z = -100;
		this.camera.position.x = 0;
		this.cameraTarget = new THREE.Vector3(0, 0,0);
		this.camera.lookAt(this.cameraTarget);
	
    const ambientLight = new THREE.AmbientLight(0xffffff, 2.485434543257532104);
    this.scene.add(ambientLight);	
	this.renderer = new THREE.WebGLRenderer({ antialias: true });
	this.renderer.setClearColor(0x1e3278, 1);
	// controls = new OrbitControls(this.camera, this.renderer.domElement);
    //controls.enabled = false
	
	
	var m = new THREE.MeshPhongMaterial({ color: 0xff0000})
	var g = new THREE.SphereGeometry(10, 6, 6);
	var mesh = new THREE.Mesh(g, m);
	mesh.position.x=0
	mesh.position.y=0
	mesh.position.z=0
	this.scene.add(mesh);
    this.renderer.setSize(200,200);
        let mapDiv = ui.div([], 'd4-viewer-host');
        this.mapDiv = mapDiv;
				mapDiv.appendChild(this.renderer.domElement);
console.error(this.renderer);
        this.root.appendChild(mapDiv);
       // mapDiv.appendChild(ms);
        //mapDiv.appendChild(this.canvas);
		this.render();
      }  
    
	render() {
		this.renderer.render(this.scene, this.camera);
	}
      
    
    
    }