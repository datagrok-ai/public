import * as THREE from 'three'
//import { OrbitControls } from "https://threejs.org/examples/jsm/controls/OrbitControls.js";
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
export class ScatterPlot3Dviewer extends DG.JsViewer {
	constructor() {
		super();
		console.log("HREEE ", THREE)

		








		this.initLayout();

		//this.onTableAttached();
	}

	getMyLook() {
		this.look.TETRAHEDRON = 'tetrahedron';
		this.look.OCTAHEDRON = 'octahedron';
		this.look.CYLINDER = 'cylinder';
		this.look.DODECAHEDRON = 'dodecahedron';
		this.look.BOX = 'box';
		this.look.SPHERE = 'sphere';
		this.look.MarkerTypes = [
			this.look.TETRAHEDRON, this.look.OCTAHEDRON, this.look.CYLINDER, 
			this.look.DODECAHEDRON, this.look.BOX, this.look.SPHERE
		];
		this.look.backColor = 'red';
		this.look.filteredRowsColor = 'green';
		this.look.filteredOutRowsColor = 'red';
		this.look.selectedRowsColor = 'green';
		this.look.missingValueColor = 'red';
		this.look.axisLineColor = 'green';
		this.look.axisTextColor = 'red';
		this.look.gridLineColor = 'green';
		/*
List<int> linearColorScheme = Color.schemeBlueWhiteRed;
  List<int> categoricalColorScheme = Color.defaultList;

		*/
		this.look.dynamicCameraMovement = false;
		this.look.showVerticalGridLines = true;
		this.look.showHorizontalGridLines = true;
		this.look.showXAxis = true;
		this.look.showYAxis = true;
		this.look.showZAxis = true;
		this.look.showFilteredOutPoints = false;

		/// Highlight 'mouse-over' rows (such as the ones that fall into a histogram bin that
  		/// the mouse is currently hovering over).
  		this.look.showMouseOverRowGroup = true;
		this.look.markerType = ScatterPlot3dMarkers.OCTAHEDRON;
		this.look.markerTypeChoices = this.look.markersTypes;
		this.look.markerRandomRotation = false;
	}

	rendererResize(size) {
		this.renderer.setSize(size.width, size.height);
	}

	initLayout() {
		this.camera = new THREE.PerspectiveCamera(45, 1, 2, 100000);
		this.scene = new THREE.Scene();
		this.camera.position.z = -100;
		this.camera.position.x = 0;
		this.cameraTarget = new THREE.Vector3(0, 0, 0);
		this.camera.lookAt(this.cameraTarget);

		const ambientLight = new THREE.AmbientLight(0xffffff, 2.485434543257532104);
		this.scene.add(ambientLight);
		this.renderer = new THREE.WebGLRenderer({ antialias: true });
		this.renderer.setClearColor(0x1e3278, 1);


		var m = new THREE.MeshPhongMaterial({ color: 0xff0000 })
		var g = new THREE.SphereGeometry(2, 6, 6);
		var mesh = new THREE.Mesh(g, m);
		mesh.position.x = 0
		mesh.position.y = 0
		mesh.position.z = 0
		this.scene.add(mesh);
		this.renderer.setSize(400, 400);
		let mapDiv = ui.div([], 'd4-viewer-host');
		this.mapDiv = mapDiv;
		mapDiv.appendChild(this.renderer.domElement);
		this.controls = new OrbitControls(this.camera, this.renderer.domElement);
		this.controls.enabled = true

		console.error(this.renderer);
		this.root.appendChild(mapDiv);
		// mapDiv.appendChild(ms);
		//mapDiv.appendChild(this.canvas);
		this.render();
	}

	// normalization of array, result: 
	// 1. array with values from 0 to 1
	// 2. scale factor (how much original array bigger than nomalized)
	// 3. center of original array
	normalize(ar) {
		var { center, scale, min, max } = this.getNormalizeInfo(ar);
		var rez = ar.map(e => (e - min) / scale)
		return rez;
	}

	// get info about array:
	// center, scale, min, max
	getNormalizeInfo(ar) {
		var min = 1000 * 1000 * 1000;
		var max = -1000 * 1000 * 1000;
		for (var i = 0; i < ar.length; i++) {
			if (ar[i] > max) max = ar[i];
			if (ar[i] < min) min = ar[i];
		};
		var dx = max - min;
		return {
			center: (max - min) / 2,
			scale: dx,
			min: min,
			max: max
		}
	}

	createMarker(type, coords) {
		var m = new THREE.MeshPhongMaterial({ color: 0xff0000 })
		var g = new THREE.SphereGeometry(2, 6, 6);
		var mesh = new THREE.Mesh(g, m);
		mesh.position.x = coords[0];
		mesh.position.y = coords[1];
		mesh.position.z = coords[2];
		return mesh
	}

	onTableAttached() {
		var df = this.dataFrame;
		// debugger 


		let numericalColumns = Array.from(this.dataFrame.columns.numerical);
		this.xColumnName = numericalColumns[0].name;
		this.yColumnName = numericalColumns[1].name;
		this.zColumnName = numericalColumns[2].name;
		this.zColumnName = 'AGE'
		this.xColumnName = 'WEIGHT';
		this.yColumnName = 'HEIGHT';
		this.rawX = this.dataFrame.getCol(this.xColumnName).getRawData();
		this.rawY = this.dataFrame.getCol(this.yColumnName).getRawData();
		this.rawZ = this.dataFrame.getCol(this.zColumnName).getRawData();
		//	console.error('raw X ', this.rawX)
		//	console.error('raw Y ', this.rawY)
		//	console.error('raw Z ', this.rawZ)
		this.placeMarkers();
		var { centerX } = this.getNormalizeInfo(this.rawX)
		var { centerY } = this.getNormalizeInfo(this.rawY)
		var { centerZ } = this.getNormalizeInfo(this.rawZ)
		this.camera.lookAt(new THREE.Vector3(centerX, centerY, centerZ))

		this.render()
	}

	placeMarkers() {
		for (var i = 0; i < this.rawX.length; i++) {
			var marker = this.createMarker('circle',
				[this.rawX[i], this.rawY[i], this.rawZ[i]]);
			this.scene.add(marker);
		}
	}

	render() {
		if (this.controls) this.controls.update();
		this.renderer.render(this.scene, this.camera);
		requestAnimationFrame(this.render.bind(this))
	}


}