/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {openPymol} from './patinae-view';
export * from './package.g';

export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.fileViewer({fileViewer: 'pse,prs,pml'})
  static previewPymol(file: DG.FileInfo): DG.View {
    return openPymol(file).view;
  }
}
