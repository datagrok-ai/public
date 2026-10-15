/* Geometry claims: where an element stands against another, how tall it is, whether it is cut off,
   and the layout audit over its subtree (`src/runtime/layout.ts`). Each claim polls the box of the one
   visible match of each phrase, allows 1 px of rounding and fails with both boxes in one sentence.
   Steps only — nouns live in kinds.ts. */
import {Page} from '@playwright/test';
import {Then} from '../../src/registry.js';
import type {ElementRef} from '../../src/runtime/args.js';
import {expectHeight, expectNoLayoutProblems, expectNotClipped, expectPlaced} from '../../src/runtime/layout.js';

export const leftOf = Then('{element} should be left of {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'left of', b), {description: 'its right edge at or before the left edge of the other'});
export const rightOf = Then('{element} should be right of {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'right of', b), {description: 'its left edge at or after the right edge of the other'});
export const above = Then('{element} should be above {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'above', b), {description: 'its bottom at or above the top of the other'});
export const below = Then('{element} should be below {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'below', b), {description: 'its top at or below the bottom of the other'});
export const inside = Then('{element} should be inside {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'inside', b), {description: 'every edge within the box of the other'});
export const clearOf = Then('{element} should not overlap {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'clear of', b), {description: 'the two boxes share at most 1 px in either direction'});
export const leftAligned = Then('{element} should be left-aligned with {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'left-aligned with', b), {description: 'the left edges of the two boxes within 1 px'});
export const topAligned = Then('{element} should be top-aligned with {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'top-aligned with', b), {description: 'the top edges of the two boxes within 1 px'});
export const centeredWith = Then('{element} should be vertically centered with {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'vertically centered with', b), {description: 'the middles of the two boxes within 1 px vertically — a label beside its editor'});
export const asTallAs = Then('{element} should be as tall as {element}', (page: Page, a: ElementRef, b: ElementRef) =>
  expectPlaced(page, a, 'as tall as', b), {description: 'the heights of the two boxes within 1 px — an input beside the Dart one it matches'});
export const tallPx = Then('{element} should be {int} px tall', (page: Page, target: ElementRef, height: number) =>
  expectHeight(page, target, height), {description: 'its box height within 1 px'});
export const notClipped = Then('{element} should not be clipped', (page: Page, target: ElementRef) => expectNotClipped(page, target),
  {description: 'whole inside every ancestor that clips it, and its own text not cut or ellipsized'});
export const noLayoutProblems = Then('{element} should have no layout problems', (page: Page, target: ElementRef) =>
  expectNoLayoutProblems(page, target),
{description: 'the layout audit over its u2 controls finds nothing: no control without size, out of reach or lying on another, ' +
  'no text cut or ellipsized, no content spilling out, labels on their editors, form editors in one column; fails on an element with no u2 control'});
