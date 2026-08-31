/* The All inputs bench as the full catalogue (packages/U2Demo/src/convergence.ts): every row
   exists, each either renders a real Dart editor beside its u2 cell or stands as an in-place
   `u2demo-ab-note` naming the platform gap, and building the page raises no console error.
   Read-only on purpose: nothing here opens a popup or edits a cell (a failing assertion under an
   open overlay OOMs the runner). */
import {ok, shot, consoleErrors} from '../local.mjs';
import {openDemoPage} from '../lib.mjs';

const GRID = '.u2demo-ab';

/** The nine Ask-H rows plus the 2.7 dataframe-picker family. */
const NEW_ROWS = ['query', 'readme', 'script', 'attachments', 'outputDir', 'privateKey',
  'owner', 'reviewers', 'param', 'table', 'keyColumn', 'features', 'mapping', 'tables'];

let errorsBefore = 0;

export async function fixture(page) {
  errorsBefore = consoleErrors.length;
  await openDemoPage(page, 'all-inputs');
  await page.waitForSelector(`${GRID} [data-row="query"]`, {timeout: 30000});
}

async function checkRowCount(page) {
  const names = await page.evaluate((grid) =>
    [...new Set([...document.querySelectorAll(`${grid} [data-row]`)].map((el) => el.dataset.row))],
  GRID);
  ok('u2demo-inputs/1 the catalogue holds at least 30 distinct rows',
    names.length >= 30, `distinct data-row names=${names.length}`);
  if (names.length < 30)
    await shot(page, 'u2demo-inputs-1-row-count');
}

async function checkNewRowsPresent(page) {
  const missing = await page.evaluate(({grid, names}) =>
    names.filter((n) => document.querySelector(`${grid} [data-row="${n}"]`) == null),
  {grid: GRID, names: NEW_ROWS});
  ok('u2demo-inputs/2 every WO-2 row name is present in the grid',
    missing.length === 0, `missing=${JSON.stringify(missing)}`);
  if (missing.length > 0)
    await shot(page, 'u2demo-inputs-2-missing-rows');
}

async function checkRowsRenderOrCarryNote(page) {
  const bad = await page.evaluate(({grid, names}) => names.map((n) => {
    const els = [...document.querySelectorAll(`${grid} [data-row="${n}"]`)];
    const isNote = els.some((el) => el.classList.contains('u2demo-ab-note'));
    // a live row is 3 cells ('…' button, Dart editor, u2 editor); a platform gap is one note
    const dartCell = els.some((el) =>
      [...el.classList].some((c) => c.startsWith('ui-input') || c === 'ui-files-input'));
    const pass = isNote ? els.length === 1 : els.length === 3 && dartCell;
    return pass ? null : `${n}: cells=${els.length} note=${isNote} dart=${dartCell}`;
  }).filter((x) => x != null), {grid: GRID, names: NEW_ROWS});
  ok('u2demo-inputs/3 each new row renders both editors or stands as a note',
    bad.length === 0, bad.join(' · '));
  if (bad.length > 0)
    await shot(page, 'u2demo-inputs-3-broken-rows');
}

async function checkConsoleClean(page) {
  const raised = consoleErrors.slice(errorsBefore);
  ok('u2demo-inputs/4 building the full catalogue raised no console error',
    raised.length === 0, raised.slice(0, 3).join(' · '));
  if (raised.length > 0)
    await shot(page, 'u2demo-inputs-4-console');
}

export const checks = [
  {id: 'u2demo-inputs/1 row count', run: checkRowCount},
  {id: 'u2demo-inputs/2 new rows present', run: checkNewRowsPresent},
  {id: 'u2demo-inputs/3 rows render or carry a note', run: checkRowsRenderOrCarryNote},
  {id: 'u2demo-inputs/4 console clean', run: checkConsoleClean},
];
