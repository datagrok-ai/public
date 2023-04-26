import {FirstViewer} from './first-viewer';
import {SecondViewer} from './second-viewer';

//tags: viewer
//output: viewer v
//meta.icon: images/1.png
//meta.toolbox: true
export function _FirstViewer() {
  return new FirstViewer();
}

//tags: viewer
//output: viewer v
export function _SecondViewer() {
  return new SecondViewer();
}

