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
    tags: ['viewer'],
    name: 'TestViewerForProperties',
    description: 'Viewer to test properties and others',
    outputs: [{type: 'viewer', name: 'result'}],
  })
  static testViewerForProperties() : DG.Viewer {
    return new TestViewerForProperties();
  }


  // -- Filters --

  @grok.decorators.func({
    tags: ['filter'],
    description: 'Test custom filter',
    outputs: [{type: 'filter', name: 'result'}],
  })
  static testCustomFilter(): DG.Filter {
    return new TestCustomFilter();
  }
}
