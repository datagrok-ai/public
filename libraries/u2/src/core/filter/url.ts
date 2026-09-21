import {computed} from '../signals.js';
import type {ReadonlySignal} from '../signals.js';

/** `base?q=<query>`, or `base` alone for an empty query; a signal in gives a signal out. */
export function queryPath(base: string, query: string): string;
export function queryPath(base: string, query: ReadonlySignal<string>): ReadonlySignal<string>;
export function queryPath(base: string, query: string | ReadonlySignal<string>): string | ReadonlySignal<string> {
  if (typeof query !== 'string')
    return computed(() => queryPath(base, query.value) as string);
  return query === '' ? base : `${base}?q=${encodeURIComponent(query)}`;
}
