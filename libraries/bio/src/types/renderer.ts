import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';

export interface IRenderer {
  get onRendered(): Observable<void>;

  invalidate(caller?: string): void;

  /** @param timeout  Default 5000 ms */
  awaitRendered(timeout?: number): Promise<void>;
}

export function isRenderer(value: IRenderer): value is IRenderer {
  return value && value.onRendered !== undefined &&
    value.invalidate !== undefined && value.awaitRendered !== undefined;
}

export function asRenderer(o: any): IRenderer | null {
  return isRenderer(o) ? o : null;
}
