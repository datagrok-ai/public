import * as THREE from 'three'
import { Vector3 } from 'three';
//import { OrbitControls } from "https://threejs.org/examples/jsm/controls/OrbitControls.js";
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls.js';

export class ScatterPlot3Dviewer extends DG.JsViewer {
	constructor() {

		super();
		console.log("TTHREEE ", THREE.REVISION);
		function hName() {
	//		console.log('hName hit!')
		}
		this.new = true;
 
		window.onHitHandlerName = 'hName';
		this.onHitHandlerName = 'hName'
 
		window.hName = hName;
		this.hName = hName
  
		//this.onHitHandlerName = onHitHandlerName;
 
		this.currentRow = -1;
		this.mouseOverRow = -1;
		this.showTooltip = true;
     
		if (this.new) {
			this.getMyLook();         
			this.initLayout();
 
			this.init3D();
//return
			console.log(this.renderer)
			if (this.isTableAttached) {
				console.error(this.dataFrame)
				this.onTableAttached()
			this.updateAllScene(this.look, this.rawX, this.rawY, this.rawZ, this.dataFrame.filter,
				this.rawZ, this.rawZ)
			}	
		} else {
			this.initLayout2()


		}

		//this.onTableAttached(); 
	} // ctor 

	rendererResize(size) {
		//	this.renderer.setSize(size.width, size.height);
	}
 
	initLayout21() {
		this.camera = new THREE.PerspectiveCamera(45, 1, .3, 100000);
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
		//mapDiv.add(ui.h1('Hello World'))
		mapDiv.appendChild(this.renderer.domElement);
		
		//this.controls.enabled = true


 
		console.error(this.renderer);
		this.root.appendChild(mapDiv);
		//this.controls = new TrackballControls(this.camera, this.renderer.domElement		);
		this.controls = new OrbitControls(this.camera, this.renderer.domElement);
		this.controls.enabled = true
		
		// mapDiv.appendChild(ms);
		//mapDiv.appendChild(this.canvas);
		this.render();
	} // initLayout2
       
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
		var m = new THREE.MeshPhongMaterial({ color: 0xff00ff })
		var g = new THREE.SphereGeometry(2, 6, 6);
		var mesh = new THREE.Mesh(g, m);
		mesh.position.x = coords[0];
		mesh.position.y = coords[1];
		mesh.position.z = coords[2];
		return mesh
	}

 
	// -------------------------------------------------------------------------------
	onTableAttached() {
		if (!this.scene) {
			this.isTableAttached = true;
			return 0;
		}
		var df = this.dataFrame;
   
//return 0            
       
		let numericalColumns = Array.from(this.dataFrame.columns.numerical);
		this.xColumnName = numericalColumns[0].name;
		this.yColumnName = numericalColumns[1].name;
		this.zColumnName = numericalColumns[2].name;
	//	this.zColumnName = 'AGE'
	//	this.xColumnName = 'WEIGHT';
	//	this.yColumnName = 'HEIGHT';
		this.rawX = this.dataFrame.getCol(this.xColumnName).getRawData();
		this.rawY = this.dataFrame.getCol(this.yColumnName).getRawData();
		this.rawZ = this.dataFrame.getCol(this.zColumnName).getRawData();
//			console.error('raw X ', this.rawX)
	//		console.error('raw Y ', this.rawY)
	//		console.error('raw Z ', this.rawZ)
//		this.placeMarkers();
		var { centerX } = this.getNormalizeInfo(this.rawX)
		var { centerY } = this.getNormalizeInfo(this.rawY)
		var { centerZ } = this.getNormalizeInfo(this.rawZ)
		this.camera.lookAt(new THREE.Vector3(centerX, centerY, centerZ))

		if (this.new) {
	//		this.placeMarkers()
	//		console.log(this.scene)
	//		this.scene.background = new THREE.Color(.2,0,0)
	//		this.render()
//return 0
			console.log('table attach if this new')
			this.clean()
		//	this.render()
	//		return 0

			this.updateAllScene(this.look, this.rawX, this.rawY, this.rawZ, this.dataFrame.filter,
				this.rawZ, this.rawZ)

	
		} else {
			this.placeMarkers()
		}
		this.render()
	} // table attached

	placeMarkers() {
		for (var i = 0; i < this.rawX.length; i++) {
			var marker = this.createMarker('circle',
				[this.rawX[i], this.rawY[i], this.rawZ[i]]);
			this.scene.add(marker);
		}
	}

	renderttt() {
		if (this.controls) this.controls.update();
		this.renderer.render(this.scene, this.camera);
		requestAnimationFrame(this.render.bind(this))
	}








	
	initLayout() {		
		this.container = ui.div([], 'd4-viewer-host');
		this.container.style.width = '100%';
		this.container.style.height = '100%';
		console.log('cont ', this.container);
	}
	
	initTHREE() {
		this.camera = new THREE.PerspectiveCamera(45, 1, .3, 100000);
		this.scene = new THREE.Scene();
		this.camera.position.z = -10;
		this.camera.position.x = 0;
		this.cameraTarget = new THREE.Vector3(0, 0, 0);
		this.camera.lookAt(this.cameraTarget);

		const ambientLight = new THREE.AmbientLight(0xffffff, 2.485434543257532104);
		this.scene.add(ambientLight);
		this.renderer = new THREE.WebGLRenderer({ antialias: true });
		this.container.appendChild(this.renderer.domElement);
		this.root.appendChild(this.container)
		this.renderer.setClearColor(0x1e3278, 1);
	}

	resetCamera2() {
		this.camera = new THREE.PerspectiveCamera(45, 1, .3, 100000);

		this.camera.position.z = -100;
		this.camera.position.x = 0;
		this.cameraTarget = new THREE.Vector3(0, 0, 0);
		this.camera.lookAt(this.cameraTarget);
	}
            
	init3D() {
		const width = this.container.clientWidth;
		const height = this.container.clientHeight;

		this.materialList = [];
		this.geometryList = [];
		this.geometrySize = new THREE.Vector3();
		this.mouse = new THREE.Vector2();
		this.highlightBoxScale = 1.5;

		//this.initTHREE();		return 
		// camera
		this.camera = new THREE.PerspectiveCamera();
		this.resetCamera();
		this.clean();
 
		
		const ambientLight = new THREE.AmbientLight(0xffffff, 2.485434543257532104);
		this.scene.add(ambientLight);

		// hitting render target
		this.hittingRenderTarget = new THREE.WebGLRenderTarget(width, height);
		this.hittingRenderTarget.texture.generateMipmaps = false;
		this.hittingRenderTarget.texture.minFilter = THREE.NearestFilter;

		var highlightBox = function (color, transparent, opacity) {
			return new THREE.Mesh(
				new THREE.BoxGeometry(1, 1, 1),
				new THREE.MeshLambertMaterial({
					emissive: color,
					transparent: transparent,
					opacity: opacity,
					side: THREE.FrontSide
				})
			);
		};
 
		// highlight boxes
		this.mouseOverBox = highlightBox(0xFFFF00, true, 0.5);
		this.currentRowBox = highlightBox(0x000000, true, 0.75);

		// renderer
		this.renderer = new THREE.WebGLRenderer({
	//		antialias: true,
		//	alpha: true
		});

//		if (this.renderer.extensions.get('ANGLE_instanced_arrays') === null) {
//			document.getElementById("notSupported").style.display = "";
//			return;
//		}


//		this.renderer.setPixelRatio(window.devicePixelRatio);
//		this.renderer.setSize(width, height);
		this.renderer.domElement.style.width = '100%';
		this.renderer.domElement.style.height = '100%';
	//	this.renderer.domElement.width = '22';
	//	this.renderer.domElement.height = '22';

		this.container.appendChild(this.renderer.domElement);
 

//		let mapDiv = ui.div([], 'd4-viewer-host');
//		this.mapDiv = mapDiv;
	//	mapDiv.add(ui.h1('Hello World'))
//		mapDiv.appendChild(this.renderer.domElement);


		this.root.appendChild(this.container)
	//	return
	//	if (this.renderer.extensions.get('instanced_arrays') === null) {
	//		throw 'ANGLE_instanced_arrays not supported';
	//	}

		// controls
		this.controls = new TrackballControls(
			this.camera, this.renderer.domElement
		);
		this.controls.dynamicDampingFactor = 0.00001;
		this.timerAutoRotation = null;

		// shaders
		this.initShaders();

		this.isEventsLinked = false;
	} // init3d


	getMyLook() {
		this.look = {}
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
		this.look.backColor = 'blue';
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
		//this.look.markerType = ScatterPlot3dMarkers.OCTAHEDRON;
		this.look.markerType = this.look.BOX;
		this.look.markerTypeChoices = this.look.markersTypes;
		this.look.markerRandomRotation = false;


		this.look.markerMinSize = 0.001;
	} // look


	updateAllScene(look, x, y, z, filter, colors, sizes) {
		this.look = look;
		this.initMesh(x, y, z, filter, colors, sizes);

		if (this.look.dynamicCameraMovement)
			this.enableAutoRotation();
		else
			this.disableAutoRotation();

		if (!this.isEventsLinked) {
			this.renderer.domElement.addEventListener('mousemove', this.onMouseMove.bind(this));
			this.renderer.domElement.addEventListener('mousedown', this.onMouseDown.bind(this));
			this.renderer.domElement.addEventListener('wheel', this.onMouseMove.bind(this));
			this.renderer.domElement.addEventListener('dblclick', this.onDoubleClick.bind(this));
			//onResize(this.container, this.onWindowResize.bind(this)); 
			requestAnimationFrame(this.render.bind(this));
			this.enableAutoRotation();
			this.isEventsLinked = true;
		}
	} 

 
	update(colors, filter) {
		let markerColor = new THREE.Color();
		for (let n = 0; n < colors.length; n++) {
			markerColor.setHex(colors[n]);
			this.markerColors.setXYZ(n, markerColor.r, markerColor.g, markerColor.b);
			this.markerAlphas.setX(n, ((colors[n] >> 24) & 255) / 255.0);
			this.filter.setX(n, filter[n]);
		}

		this.igeo.attributes.color.needsUpdate = true;
		this.igeo.attributes.alpha.needsUpdate = true;
		this.igeo.attributes.filter.needsUpdate = true;
		this.hit();
		this.render();
	}



	hitRows(currentRow, mouseOverRow) {
		this.currentRow = currentRow;
		this.mouseOverRow = mouseOverRow;
		this.hit();
		this.render();
	}



	resetView() {
		this.resetCamera();
		this.controls.reset();
		this.hit();
		this.render();
	}


	enableAutoRotation() {
		return 0;
		this.controls.staticMoving = false;

		// Simulate touch
		this.controls.handleEvent(new MouseEvent('mousedown', { clientX: 0 }));
		this.controls.handleEvent(new MouseEvent('mousemove', { clientX: 1 }));
		this.controls.handleEvent(new MouseEvent('mouseup', { clientX: 2 }));

		if (this.timerAutoRotation === null)
			this.timerAutoRotation = setInterval(this.render.bind(this), 10);
	}


	disableAutoRotation() {
		this.controls.staticMoving = true;
		clearInterval(this.timerAutoRotation);
		this.timerAutoRotation = null;
	}


	resetCamera() {
		this.camera.fov = 45;
		this.camera.aspect =1// this.container.clientWidth / this.container.clientHeight;
		this.camera.near = 1;
		this.camera.far = 100;
		this.camera.position.z = 4.0;
	//	this.camera.updateProjectionMatrix();
	}


	initMesh(x, y, z, filter, color, size) {
		this.clean();
		console.log('init mesh', this.look.markerType)
		let geo;
		if (this.look.markerType === 'tetrahedron')
			geo = new THREE.TetrahedronGeometry(1.0);
		else if (this.look.markerType === 'octahedron')
			geo = new THREE.OctahedronGeometry(1.0);
		else if (this.look.markerType === 'cylinder')
			geo = new THREE.CylinderGeometry(1.0, 1.0, 1.0, 5, 1);
		else if (this.look.markerType === 'dodecahedron')
			geo = new THREE.DodecahedronGeometry(1.0);
		else if (this.look.markerType === 'box')
			geo = new THREE.BoxGeometry(1.0, 1.0, 1.0);
		else if (this.look.markerType === 'sphere')
			geo = new THREE.SphereGeometry(1.0, 8, 8);
		else
			geo = new THREE.TetrahedronGeometry(1.0);
		geo.computeBoundingBox();
		geo.boundingBox.getSize(this.geometrySize);
		this.geometryList.push(geo);

		this.makeInstanced(geo, x, y, z, filter, color, size);
		this.hit();
		this.render();
	}

	initMesh2() {

	}


	saveMousePosition(e) {
		this.mouse.x = e.clientX;
		this.mouse.y = e.clientY;
	}


	onMouseMove(e) {
		if (this.showTooltip === false)
		debugger
			window[this.onHitHandlerName].apply(window, [-1, this.currentRow, 0, 0, new Object()]);
		this.saveMousePosition(e);
		const element = this.hitTest(e.layerX, e.layerY);
		this.mouseOverRow = element.id;
		this.hit();
		if (!this.look.dynamicCameraMovement)
			this.render();
	}



	onMouseDown(e) {
		if (!this.look.dynamicCameraMovement)
			this.disableAutoRotation();
		handleMouseMove(this.onMouseMove.bind(this), this.onMouseUp.bind(this));
		this.showTooltip = false;
		this.saveMousePosition(e);
		const element = this.hitTest(e.layerX, e.layerY);
		if (element.object && !(e.ctrlKey || e.metaKey || e.shiftKey))
			this.currentRow = element.id;

		let options = {
			'ctrlKey': e.ctrlKey,
			'metaKey': e.metaKey,
			'shiftKey': e.shiftKey
		};

		this.hit(options);
		this.render();
	}


	onMouseUp(e) {
		this.showTooltip = true;
	}


	onDoubleClick(e) {
		const element = this.hitTest(e.layerX, e.layerY);
		if (element.object)
			this.onMouseDown(e);
		else
			this.resetView();
	}




	onWindowResize(e) {
		const width = this.container.clientWidth;
		const height = this.container.clientHeight;
		this.camera.aspect = width / height;
		this.camera.updateProjectionMatrix();

		this.renderer.setSize(width, height);
		this.hittingRenderTarget.setSize(width, height);
		this.renderer.render(this.scene, this.camera);
	}




	clean() {
	/*	
		this.scene = new THREE.Scene();
		this.scene.background = new THREE.Color(this.look.backColor);
		this.scene.background = new THREE.Color(0,0,1);
		return 0;
*/
		THREE.Cache.clear();

		this.materialList.forEach(function (m) {
			m.dispose();
		});

		this.geometryList.forEach(function (g) {
			g.dispose();
		});

		this.scene = new THREE.Scene();
		this.scene.background = new THREE.Color(this.look.backColor);
		this.scene.background = new THREE.Color(0, 1, 0);

		this.scene.add(this.camera);
	//	this.scene.add(this.mouseOverBox);
	//	this.scene.add(this.currentRowBox);

		this.hittingScene = new THREE.Scene();
		this.hittingData = {};
		this.materialList = [];
		this.geometryList = [];
	}





	makeInstanced(geo, x, y, z, filter, colors, sizes) {
	//	return 0;
		console.log('make instanced ', geo)
		// material
		const vert = document.getElementById('vertInstanced').textContent;
		const frag = document.getElementById('fragInstanced').textContent;

		const material = new THREE.RawShaderMaterial({
			vertexShader: vert,
			fragmentShader: frag,
			transparent: true
		});

		this.materialList.push(material);

		const hittingMaterial = new THREE.RawShaderMaterial({
			vertexShader: "#define HITTING\n" + vert,
			fragmentShader: "#define HITTING\n" + frag
		});

		this.materialList.push(hittingMaterial);

		// geometry 
		console.log('geo ', geo);   
		//let bg2 =  new THREE.BufferGeometry(geo)
		//console.log('bg2 ', bg2)
		//let bgeo =bg2.fromGeometry(geo);
		//let bgeo =bg2.fromGeometry(geo);
		//let bgeo =bg2
		var bgeo = geo.clone();
		this.geometryList.push(bgeo);

		this.igeo = new THREE.InstancedBufferGeometry();
		this.geometryList.push(this.igeo);

		const vertices = bgeo.attributes.position.clone();
		this.igeo.addAttribute('position', vertices);

		var numObjects = x.length;
		numObjects = 100

		let mcol0 = new THREE.InstancedBufferAttribute(new Float32Array(numObjects * 3), 3, false);
		let mcol1 = new THREE.InstancedBufferAttribute(new Float32Array(numObjects * 3), 3, false);
		let mcol2 = new THREE.InstancedBufferAttribute(new Float32Array(numObjects * 3), 3, false);
		let mcol3 = new THREE.InstancedBufferAttribute(new Float32Array(numObjects * 3), 3, false);

		const scaleFactor = 1.0 / (10.0 * Math.log(numObjects));

		let position = new THREE.Vector3();
		let rotation = new THREE.Euler();
		let quaternion = new THREE.Quaternion();
		let scale = new THREE.Vector3();

		let matrix = new THREE.Matrix4();
		let me = matrix.elements;

		//for (let n = 0; n < numObjects; n++) {
		for (let n = 0; n < numObjects; n++) {
			let size = (sizes === null) ? scaleFactor : sizes[n] * scaleFactor;
			size = (size < this.markerMinSize) ? this.markerMinSize : size;
			scale.setScalar(size);

			position.x = x[n] - 0.5;
			position.y = y[n] - 0.5;
			position.z = z[n] - 0.5;

			if (this.look.markerRandomRotation) {
				rotation.x = Math.random() * 2 * Math.PI;
				rotation.y = Math.random() * 2 * Math.PI;
				rotation.z = Math.random() * 2 * Math.PI;
				quaternion.setFromEuler(rotation, false);
			}
   
			matrix.compose(position, quaternion, scale);

			let object = new THREE.Object3D();
			object.applyMatrix4(matrix); 
			this.hittingData[n + 1] = object;

			mcol0.setXYZ(n, me[0], me[1], me[2]);
			mcol1.setXYZ(n, me[4], me[5], me[6]);
			mcol2.setXYZ(n, me[8], me[9], me[10]);
			mcol3.setXYZ(n, me[12], me[13], me[14]);
		}
/*
		this.igeo.addAttribute('mcol0', mcol0);
		this.igeo.addAttribute('mcol1', mcol1);
		this.igeo.addAttribute('mcol2', mcol2);
		this.igeo.addAttribute('mcol3', mcol3);
*/
this.igeo.setAttribute('mcol0', mcol0);
this.igeo.setAttribute('mcol1', mcol1);
this.igeo.setAttribute('mcol2', mcol2);
this.igeo.setAttribute('mcol3', mcol3);
		this.markerColors = new THREE.InstancedBufferAttribute(new Float32Array(numObjects * 3), 3, false);
		this.markerAlphas = new THREE.InstancedBufferAttribute(new Float32Array(numObjects), 1, false);
		let markerColor = new THREE.Color();

		for (let n = 0; n < numObjects; n++) {
			markerColor.setHex(colors[n]);
			this.markerAlphas.setX(n, ((colors[n] >> 24) & 255) / 255.0);
			this.markerColors.setXYZ(n, markerColor.r, markerColor.g, markerColor.b);
		}

		this.igeo.setAttribute('color', this.markerColors);
		this.igeo.setAttribute('alpha', this.markerAlphas);

		let col = new THREE.Color();
		let hittingColors = new THREE.InstancedBufferAttribute(new Float32Array(numObjects * 3), 3, false);

		for (let n = 0; n < numObjects; n++) {
			col.setHex(n + 1);
			hittingColors.setXYZ(n, col.r, col.g, col.b);
		}

		this.igeo.setAttribute('hittingColor', hittingColors);

		// filter
		this.filter = new THREE.InstancedBufferAttribute(filter, 1, false);
		this.igeo.setAttribute('filter', this.filter);

		// mesh
		var m = new THREE.MeshPhongMaterial({ color: 0x0000ff })

		this.mesh = new THREE.Mesh(this.igeo, m);
		this.mesh.name = 'name'
		this.mesh.scale.x =100
		this.mesh.scale.y =100
		console.log('make instanced ', geo)

		console.log('mesh ', this.mesh);
		this.scene.add(this.mesh);




		this.hittingMesh = new THREE.Mesh(this.igeo, hittingMaterial);
		this.hittingScene.add(this.hittingMesh);





	var m = new THREE.MeshPhongMaterial({ color: 0xff0000 })
	var g = new THREE.SphereGeometry(1, 6, 6);
	var mesh2 = new THREE.Mesh(g, m);
	mesh2.position.x = 0
	mesh2.position.y = 0
	mesh2.position.z = 0
	this.scene.add(mesh2);
	} // makeInstanced


	hitTest(x, y) {
		// create buffer for reading a single pixel
		let pixelBuffer = new Uint8Array(4);

		// read the pixel under the mouse from the texture
		this.renderer.readRenderTargetPixels(this.hittingRenderTarget,
			x, this.hittingRenderTarget.height - y, 1, 1, pixelBuffer);

		// interpret the pixel as an ID
		let id = (pixelBuffer[0] << 16) | (pixelBuffer[1] << 8) | (pixelBuffer[2]);

		return {
			id: id - 1,
			object: this.hittingData[id]
		};
	} // hitTest


	getElement(id) {
		return {
			id: id,
			object: this.hittingData[id + 1]
		};
	}


	hit(options) {
		options = ((options === undefined) || (options === null)) ? new Object() : options;

		options.showTooltip = this.showTooltip;

		// render the hitting scene off-screen
		this.mouseOverBox.visible = false;
		this.currentRowBox.visible = false;
	//	this.renderer.render(this.hittingScene, this.camera, this.hittingRenderTarget);
		this.moveHighlightBox(this.mouseOverBox, this.getElement(this.mouseOverRow), true, options);
		this.moveHighlightBox(this.currentRowBox, this.getElement(this.currentRow), false, options);
	}

 
	moveHighlightBox(highlightBox, element, notify, options) {
		if (element.object) {
			const object = element.object;

			// move the highlightBox, so that it surrounds the hit object
			if (object.position && object.rotation && object.scale) {
				highlightBox.position.copy(object.position);
				highlightBox.rotation.copy(object.rotation);
				highlightBox.scale.copy(object.scale).multiply(this.geometrySize).multiplyScalar(this.highlightBoxScale);
				highlightBox.visible = true;
			}

			if (notify)
				window[this.onHitHandlerName].apply(window, [this.mouseOverRow, this.currentRow, this.mouse.x, this.mouse.y, options]);
		}
		else {
			if (notify)
				window[this.onHitHandlerName].apply(window, [this.mouseOverRow, this.currentRow, 0, 0, options]);
		}
	}
  
      
	render() {
		if (this.controls) this.controls.update();
	//	console.log(this.controls)
		//this.camera.rotation.x = Math.PI
		if (false) {
			this.camera.position.x = 10;
			this.camera.position.y = 10; 
			this.camera.position.z = 10;
			this.camera.lookAt(new Vector3(0,0,0))
		}

		
	//	const ambientLight = new THREE.AmbientLight(0xffffff, 222.485434543257532104);
	//	this.scene.add(ambientLight); 

		  
	var m = new THREE.MeshPhongMaterial({ color: 0xff0000 })  
	var g = new THREE.SphereGeometry(1.1, 6, 6);
	var mesh2 = new THREE.Mesh(g, m);
	mesh2.position.x = 5
	mesh2.position.y = 5
	mesh2.position.z = 0
//	this.scene.add(mesh2);

	//	console.log('camera ', this.camera.position, this.camera.rotation)
	//	console.log('scene ', this.scene)
		this.renderer.render(this.scene, this.camera);
//return 0 
		requestAnimationFrame(this.render.bind(this))

	}




	initShaders() {
		const vertInstanced = `
		#define SHADER_NAME vertInstanced
  
		precision highp float;
	  
		uniform mat4 modelViewMatrix;
		uniform mat4 projectionMatrix;
	  
		attribute vec3 position;
		attribute vec3 mcol0;
		attribute vec3 mcol1;
		attribute vec3 mcol2;
		attribute vec3 mcol3;
		attribute float filter;
		attribute float alpha;
	  
		#ifdef HITTING
		attribute vec3 hittingColor;
		#else
		attribute vec3 color;
		varying vec3 vPosition;
		varying float vAlpha;
		#endif
	  
		varying vec3 vColor;
	  
		void main()	{
		  mat4 matrix = mat4(
			vec4(mcol0, 0),
			vec4(mcol1, 0),
			vec4(mcol2, 0),
			vec4(mcol3, 1)
		  );
	  
		  vec3 positionEye = (modelViewMatrix * matrix * vec4(position, 1.0)).xyz;
	  
		#ifdef HITTING
		  vColor = hittingColor;
		#else
		  vColor = color;
		  vPosition = positionEye;
		  vAlpha = alpha;
		#endif
		  gl_Position = projectionMatrix * vec4(positionEye, filter);
		  gl_Position = projectionMatrix * vec4(positionEye, 1.);
		}`;

		const fragInstanced = `
		#define SHADER_NAME fragInstanced
  
		#extension GL_OES_standard_derivatives : enable
  
		precision highp float;
	  
		varying vec3 vColor;
	  
		#ifndef HITTING
		varying vec3 vPosition;
		varying float vAlpha;
		#endif
	  
		void main()	{
		#ifdef HITTING
		  gl_FragColor = vec4(vColor, 1.0);
		#else
				  vec3 fdx = dFdx(vPosition);
				  vec3 fdy = dFdy(vPosition);
				  vec3 normal = normalize(cross(fdx, fdy));
				  float diffuse = dot(normal, vec3(0.0, 0.0, 1.0));
				  gl_FragColor = vec4(diffuse * vColor, vAlpha);
		#endif
		gl_FragColor = vec4(1., 0.,0.,1.);
		}`;

		const addScript = function (type, id, code) {
			let s = document.createElement("script");
			s.id = id;
			s.type = type;
			s.innerHTML = code;
			document.head.appendChild(s);
		};

		if (!document.getElementById("vertMerged")) {
			addScript("x-shader/x-vertex", "vertInstanced", vertInstanced);
			addScript("x-shader/x-fragment", "fragInstanced", fragInstanced);
		}

		this.vertexShader = vertInstanced;
		this.fragmentShader = fragInstanced;
	} // initShaders


}