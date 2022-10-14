import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface IViewer {
  /** JsViewer.root(): HTMLElement */
  get root(): HTMLElement;

  /** Viewer.close(): void */
  close(): void;

  /** Viewer.removeFromView(): any */
  removeFromView(): any;
}