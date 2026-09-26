import {fromEvent, Observable, timer} from 'rxjs';
import {distinctUntilChanged, filter, map, switchMap, take, timeout} from 'rxjs/operators';

/**
 * Emits (and completes) when the element returned by `get` is clicked.
 *
 * `get` is re-run on every tick rather than resolved once, so a target that is rebuilt while the
 * step is up — a ribbon item, a re-rendered toolbar — gets the listener too; binding once leaves
 * the step unfinishable on the node that was replaced.
 *
 * The timeout sits after the null filter, so it counts only while `get` returns nothing: once the
 * element exists the step waits for the click for as long as the learner needs.
 *
 * Returning an Observable lets the tutorial engine cancel the step (via `firstEvent`) when the
 * tutorial is closed.
 */
export function elementClick(get: () => HTMLElement | null, timeoutMs = 30000): Observable<Event> {
  return timer(0, 100).pipe(
    map(() => get()),
    filter((el): el is HTMLElement => el != null),
    timeout(timeoutMs),
    distinctUntilChanged(),
    switchMap((el) => fromEvent(el, 'click')),
    take(1),
  );
}
