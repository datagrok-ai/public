import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {IEditor} from '@datagrok/js-draw-lite/src/types/jsdraw2';
import {HelmType} from '@datagrok/js-draw-lite/src/types/org';

export interface IHelmWebEditor {
  get editor(): IEditor<HelmType>;
  get host(): HTMLDivElement;

  resizeEditor(width: number, height: number): void;
}
