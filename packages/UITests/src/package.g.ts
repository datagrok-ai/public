import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: TestViewerForProperties
//description: Viewer to test properties and others
//output: viewer result
//meta.role: viewer,panel
export function testViewerForProperties() : any {
  return PackageFunctions.testViewerForProperties();
}

//description: Test custom filter
//output: filter result
//meta.role: filter
export function testCustomFilter() : any {
  return PackageFunctions.testCustomFilter();
}
