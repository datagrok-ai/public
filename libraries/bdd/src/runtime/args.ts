/* The value carriers generated specs pass to step functions, and the context switch they emit
   after a step declared with `enters`. */
import type {Page} from '@playwright/test';
import {ContextEntry, DatasetEntry, lookupContext, lookupDataset} from '../registry.js';

export interface ElementRef {
  phrase: string;
}

/** An element phrase, resolved against the registry (and the page's context) when a step runs. */
export function el(phrase: string): ElementRef {
  return {phrase};
}

export function ds(name: string): DatasetEntry {
  const entry = lookupDataset(name);
  if (!entry)
    throw new Error(`unknown dataset "${name}"`);
  return entry;
}

let current: {page: Page; context: ContextEntry} | undefined;

/** From here on, phrases on this page resolve with the context's vocabulary. Bound to the page,
 * so a context never outlives its test. */
export function enter(page: Page, name: string): void {
  const context = lookupContext(name);
  if (!context)
    throw new Error(`unknown context "${name}"`);
  current = {page, context};
}

export function leave(page: Page): void {
  if (current?.page === page)
    current = undefined;
}

export function contextOf(page: Page): ContextEntry | undefined {
  return current?.page === page ? current.context : undefined;
}
