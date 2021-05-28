//const { test22func } = require("../test22");
import * from './test22'


console.log('ScatterFast3D file start ...')
class ScatterFast3dViewer extends DG.JsViewer {

  constructor() {
    console.log('ScatterFast3dViewer ctor')
    test22func();
    super();
    // properties

    this.scale = [1,1]
    this.canvasSize = [500,400]

    ui.onSizeChanged(this.root).subscribe((e) => {
      var b = e.getBoundingClientRect();
      this.updateCanvaseSize(b.width, b.height)
      if (this.gl) this.render()
    });
    this.initLayout()
  } // ctor


  onPropertyChanged(property) {
    var name = property.name;
    var val = property.get();
    if (name === 'markerSizeProp') {
      this.markerSize = parseInt(val);
      if (this.dataFrame) this.render();
    };
    //super.onPropertyChanged(property);
  }

  initLayout() {
    
    this.canvas = document.createElement('canvas');
    console.error(this.canvas)
    this.canvas.style.border = '1px solid red'
    let mapDiv = ui.div([], 'd4-viewer-host');
    this.mapDiv = mapDiv;
    this.ctx = this.canvas.getContext('2d');
    this.root.appendChild(mapDiv);
   // mapDiv.appendChild(ms);
    mapDiv.appendChild(this.canvas);
    this.canvas.style.width = '333px';
    this.canvas.style.height = '333px';
  }

  updateCanvaseSize(w, h) {
    this.canvas.width = w;
    this.canvas.height = h;
    this.canvas.style.width = '' + w + 'px';
    this.canvas.style.height = '' + h + 'px';
    this.scale[0] = this.scale[0] * w / this.canvasSize[0];
    this.scale[1] = this.scale[1] * h / this.canvasSize[1];
    this.canvasSize[0] = w;
    this.canvasSize[1] = h;
    if (this.gl) this.gl.viewport(0, 0, w, h);
  }

  onTableAttached() {
 
  } // table attached



  init() {


  } //  init



  updateStats() {
    // update fps meter
 //   this.time0 = this.time1 // save old timestamps
 //   this.time1 = Date.now()
    let fps2 = 1000 / (this.time1 - this.time0)
    this.fpsMeter.innerHTML =
      fps2.toFixed(1) + ' fps; t: ' + (this.time1-this.time0) + 
      'ms; size: ' + (this.rowCount / 1000).toFixed(0) + 'k ' +
      '<br>Scale: ' + [this.scale[0].toFixed(2), this.scale[1].toFixed(2)] +
      'pan: ' + [this.pan[0].toFixed(2), this.pan[1].toFixed(2)] + '  ' +
      '<br>size: ' + this.canvasSize + ' texHeight: ' + this.dataTextureHeight +
      '<br>Marker: ' + this.markerSize 
      //+      '<br>screenPoint: ' + this.screenPoint;
  }

  render() {


    if (this.debug) {      

      var rez = new Uint8Array(50 * 1 * 4);

      this.time1 = Date.now();

      var singleNumber = rez[0] + 256 * (rez[1] + 256 * (rez[2] + 256 * rez[3]));
      console.log('bytes ', rez[0], rez[1], rez[2], rez[3]);
      console.log('single number: ', singleNumber);
  
      this.updateStats();  
      if (this.timeControl && this.timePlayMode) {
        this.timeControlUpdate()
        window.requestAnimationFrame(this.render.bind(this))
      }
    };
  } // render
}

console.log('ffff1 e')
