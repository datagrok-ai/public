import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import '../types/jsdraw2';
import * as JSDraw2 from 'JSDraw2';

export interface IHelmWebEditor {
  get editor(): JSDraw2.Editor;
  get host(): HTMLDivElement;

  resizeEditor(width: number, height: number): void;
}
