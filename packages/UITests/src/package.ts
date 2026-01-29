/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
export * from './package.g';
import {TestViewerForProperties} from './viewers/test-viewer-for-properties';
import {TestCustomFilter} from './viewers/test-custom-filter';
export const _package = new DG.Package();


export class PackageFunctions {
  // -- Viewers --

  @grok.decorators.panel({
    name: 'TestViewerForProperties',
    description: 'Viewer to test properties and others',
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {role: 'viewer'},
  })
  static testViewerForProperties() : DG.Viewer {
    return new TestViewerForProperties();
  }


  // -- Filters --

  @grok.decorators.func({
    description: 'Test custom filter',
    outputs: [{type: 'filter', name: 'result'}],
    meta: {role: 'filter'},
  })
  static testCustomFilter(): DG.Filter {
    return new TestCustomFilter();
  }
}
