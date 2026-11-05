import * as DG from 'datagrok-api/dg';

import {defer, Observable} from 'rxjs';
import {finalize, tap} from 'rxjs/operators';

import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

/** A {@link PromiseSyncer} that knows whether a change is still on its way to the scene: a queued
 * call, or a debounced request its handler has not queued yet. What a viewer's `isRenderPending`
 * reads. */
export class PendingSyncer extends PromiseSyncer {
  private running: number = 0;
  private readonly waiting = new Set<object>();

  constructor(logger: ILogger, private readonly onIdle: () => void) {
    super(logger);
  }

  get isPending(): boolean {
    return this.running > 0 || this.waiting.size > 0;
  }

  /** {@link DG.debounce} whose request counts as pending from its first emission until the debounced
   * value reaches the handler, or until the subscription ends — a request never outlives its handler. */
  debounce<T>(source: Observable<T>, milliseconds: number): Observable<T> {
    return defer(() => {
      const request = {};
      return DG.debounce(source.pipe(tap(() => this.waiting.add(request))), milliseconds).pipe(
        tap(() => this.waiting.delete(request)),
        finalize(() => this.waiting.delete(request)));
    });
  }

  override sync(logPrefix: string, func: () => Promise<void>): void {
    this.running++;
    super.sync(logPrefix, async () => {
      try {
        await func();
      } finally {
        if (--this.running === 0)
          this.onIdle();
      }
    });
  }
}
