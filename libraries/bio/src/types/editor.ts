import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export interface IHelmWebEditor {
  get editor(): any;
  get host(): HTMLDivElement;

  createWebEditor(substructure: string): { editorDiv: HTMLDivElement, webEditor: any };
  resizeEditor(width: number, height: number): void;
}
