import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: TestViewerForProperties
//description: Viewer to test properties and others
//tags: panel, viewer
//output: viewer result
export function testViewerForProperties() : any {
  return PackageFunctions.testViewerForProperties();
}

//description: Test custom filter
//tags: filter
//output: filter result
export function testCustomFilter() : any {
  return PackageFunctions.testCustomFilter();
}
