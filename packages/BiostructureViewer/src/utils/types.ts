import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {Observable} from 'rxjs';

export interface IPdbGridCellRenderer {
  get renderComplete(): Observable<void>;
}
