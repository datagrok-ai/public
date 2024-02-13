import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IRenderer, isIRenderer} from '../types/renderer';

export interface IViewer extends IRenderer {
  /** JsViewer.root(): HTMLElement */
  get root(): HTMLElement;

  /** Viewer.close(): void */
  close(): void;

  /** Viewer.removeFromView(): any */
  removeFromView(): any;
}

export function isIViewer(value: IViewer): value is IViewer {
  return isIRenderer(value) &&
    value.root !== undefined && value.close !== undefined && value.removeFromView !== undefined;
}

export function asIViewer(o: any): IViewer | null {
  return isIViewer(o) ? o : null;
}
