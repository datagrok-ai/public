/* The pointer's approach to a target: from a point beside it, then onto it in five moves, as a person's pointer
   arrives. A tooltip that answers a moving pointer (an overlay on a disabled button, a canvas grid's cell, a scatter
   plot's point) shows; a single move does not always raise it. No steps here: the step modules import it. */
import type {Page} from '@playwright/test';

export interface Point {x: number; y: number}

export async function slidePointer(page: Page, to: Point, from: Point): Promise<void> {
  await page.mouse.move(from.x, from.y);
  await page.mouse.move(to.x, to.y, {steps: 5});
}
