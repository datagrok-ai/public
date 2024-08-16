/* Do not change these import lines to match external modules in webpack configuration */
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TestViewerForProperties} from './viewers/test-viewer-for-properties';
import {TestCustomFilter} from './viewers/test-custom-filter';
import { delay } from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();

// -- Viewers --

//name: TestViewerForProperties
//description: Viewer to test properties and others
//tags: viewer, panel
//output: viewer result
export function testViewerForProperties() {
  return new TestViewerForProperties();
}

//name: TestViewerProperties
//description: Viewer to test properties and others
//input: dataframe df
//input: string viewerType
//input: dynamic properties
//tags: viewer, panel
//output: viewer result 
export async function testViewerProperties(df : DG.DataFrame, viewerType: DG.ViewerType, properties: { [key: string]: any}) {
  var tv = DG.TableView.create(df, true);
  var viewer = tv.addViewer(viewerType);
  viewer.setOptions(properties);
  await delay(1000);
  return tv
}

// -- Filters --

//name: testCustomFilter
//description: Test custom filter
//tags: filter
//output: filter result
export function testCustomFilter(): DG.Filter {
  const flt: TestCustomFilter = new TestCustomFilter();
  return flt;
}
