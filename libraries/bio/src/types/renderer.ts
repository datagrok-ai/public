import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';

export interface IRenderer {
  get onRendered(): Observable<void>;

  invalidate(): void;

  /** @param timeout  Default 5000 ms */
  awaitRendered(timeout?: number): Promise<void>;
}
