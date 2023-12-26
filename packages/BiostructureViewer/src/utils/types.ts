import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {Observable} from 'rxjs';

import {IRenderer} from '@datagrok-libraries/bio/src/types/renderer';

export interface IPdbGridCellRenderer extends IRenderer {
  get onClicked(): Observable<void>;
}
