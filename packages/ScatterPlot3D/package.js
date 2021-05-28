
import {test22func} from './test22.js'


class ScatterFast3dPackage extends DG.Package {

  //name: ScatterFast3d
  //description: Fast Scatter Viewer 3d
  //tags: viewer app
  //output: viewer result
  scatterFast3dViewer() {
    return new ScatterFast3dViewer();
  }

  //name: tFunc
  //description: aaa
  //tags: function
  tFunc() {
    test22func();
  }
}