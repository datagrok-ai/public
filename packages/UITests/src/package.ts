/* Do not change these import lines to match external modules in webpack configuration */
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TestViewerForProperties} from './viewers/test-viewer-for-properties';
import {TestCustomFilter} from './viewers/test-custom-filter';
export const _package = new DG.Package();

// -- Viewers --

//name: TestViewerForProperties
//description: Viewer to test properties and others
//tags: viewer, panel
//output: viewer result
export function testViewerForProperties() {
  return new TestViewerForProperties();
}

// -- Filters --

//name: testCustomFilter
//description: Test custom filter
//tags: filter
//output: filter result
export function testCustomFilter(): DG.Filter {
  return new TestCustomFilter();
}
