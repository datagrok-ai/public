import {fromEvent, Observable, timer} from 'rxjs';
import {distinctUntilChanged, filter, map, switchMap, take, timeout} from 'rxjs/operators';

/**
 * Emits (and completes) when the element returned by `get` is clicked.
 * Polls for the element so it also works when it renders slightly later,
 * attaches the click listener the moment the element exists, and errors out
 * after `timeoutMs` instead of waiting forever. Returning an Observable lets the
 * tutorial engine cancel the step (via `firstEvent`) when the tutorial is closed.
 */
export function elementClick(get: () => HTMLElement | null, timeoutMs = 30000): Observable<Event> {
  return timer(0, 100).pipe(
    map(() => get()),
    filter((el): el is HTMLElement => el != null),
    timeout(timeoutMs),
    distinctUntilChanged(),
    switchMap((el) => fromEvent(el, 'click', {capture: true})),
    take(1),
  );
}
