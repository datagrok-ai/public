import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

export const MATRIX = '[name="viewer-Matrix-plot"]';
export const INNER_CELLS = `${MATRIX} canvas.d4-matrix-plot-inner-viewer`;

/**
 * Installs the onViewerRendered stamp `waitForViewerRendered` reads.
 *
 * Without it the first wait finds no stamp, so it waits for a render that already happened
 * inside the synchronous property set and burns its whole cap (0.9-1.2s, measured on every
 * matrix spec).
 */
export const stampRenders = (page: Page) => v.waitForViewerRendered(page, 'Matrix plot', 50);

interface SettleArgs {
  sel: string;
  idxs: number[];
  from: number | null;
  changeCapMs: number;
  tolerance: number;
  rounds: number;
  tickMs: number;
}

/**
 * The ink-settle loop the matrix specs each kept their own copy of, moved into the page.
 *
 * A single negative index sums every inner cell; several indices settle together, so two
 * cells cost one 300ms round rather than two. Cadence, tolerance and iteration cap are the
 * ones the Playwright-side copies used; what changes is that the reads no longer cross the
 * process boundary once per tick, and neither does the wait for the ink to move off `from`.
 */
async function settle(page: Page, args: SettleArgs): Promise<number[]> {
  return page.evaluate(async (a: SettleArgs) => {
    const one = (c: HTMLCanvasElement | undefined): number => {
      const ctx = c?.getContext('2d');
      if (!ctx) return -1;
      let data: Uint8ClampedArray;
      try { data = ctx.getImageData(0, 0, c!.width, c!.height).data; } catch (_) { return -2; }
      let n = 0;
      for (let k = 0; k < data.length; k += 16)
        if (data[k + 3] !== 0 && !(data[k] >= 250 && data[k + 1] >= 250 && data[k + 2] >= 250)) n++;
      return n;
    };
    const ink = (): number[] => {
      const cells = document.querySelectorAll(a.sel);
      return a.idxs.map((idx) => {
        if (idx >= 0) return one(cells[idx] as HTMLCanvasElement | undefined);
        let total = 0;
        for (const c of Array.from(cells) as HTMLCanvasElement[]) {
          const n = one(c);
          if (n > 0) total += n;
        }
        return total;
      });
    };
    const settled = (a2: number[], b: number[]) =>
      a2.every((x, i) => x >= 0 && Math.abs(x - b[i]) < a.tolerance);
    const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

    if (a.from !== null) {
      const deadline = Date.now() + a.changeCapMs;
      const from = a.idxs.map(() => a.from as number);
      while (Date.now() < deadline) {
        if (ink().every((x, i) => x >= 0 && Math.abs(x - from[i]) >= a.tolerance)) break;
        await sleep(100);
      }
    }
    // the criterion is unchanged — the ink now against the ink one window ago — but the window
    // slides at half its length, so an unsettled read costs half a window to notice, not a whole one
    const half = Math.round(a.tickMs / 2);
    const history = [ink()];
    for (let i = 0; i < a.rounds * 2; i++) {
      await sleep(half);
      history.push(ink());
      if (history.length > 3) history.shift();
      if (history.length === 3 && settled(history[2], history[0])) break;
    }
    return history[history.length - 1];
  }, args);
}

const cellArgs = {sel: INNER_CELLS, tolerance: 40, rounds: 10, tickMs: 300};
const matrixArgs = {sel: INNER_CELLS, idxs: [-1], tolerance: 80, rounds: 12, tickMs: 300};

export const settledCellInk = async (page: Page, idx: number): Promise<number> =>
  (await settle(page, {...cellArgs, idxs: [idx], from: null, changeCapMs: 0}))[0];

export const settledCellInks = (page: Page, idxs: number[]): Promise<number[]> =>
  settle(page, {...cellArgs, idxs, from: null, changeCapMs: 0});

export const settledCellInkAfterChange = async (
  page: Page, idx: number, from: number, capMs: number,
): Promise<number> => (await settle(page, {...cellArgs, idxs: [idx], from, changeCapMs: capMs}))[0];

export const settledMatrixInk = async (page: Page): Promise<number> =>
  (await settle(page, {...matrixArgs, from: null, changeCapMs: 0}))[0];

export const settledMatrixInkAfterChange = async (
  page: Page, from: number, capMs: number,
): Promise<number> => (await settle(page, {...matrixArgs, from, changeCapMs: capMs}))[0];

export const cellCount = (page: Page): Promise<number> =>
  page.evaluate((sel) => document.querySelectorAll(sel).length, INNER_CELLS);
