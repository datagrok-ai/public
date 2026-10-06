/* The U2Demo "Run history" pane (packages/U2Demo/src/func-history.ts): a FunctionInput on top,
   the funcForm below with the two standard history properties on — the Run button (`showRun`)
   and the history icon (`showHistory`). The local lane is serverless (LocalClient answers dapi
   with canned/empty payloads), so this checks the shape and the gating: picking a function
   through the FunctionInput popup builds the form, Run sits disabled while a required field is
   empty and arms when it is filled, and the history popup opens (empty here) and closes on Esc.
   Actually running + saving + pulling a run back is the live-stand check. */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

const PAGE = '.u2demo-func-history';
const FI = `${PAGE} [data-u2="function-input"] .u2-function-input`;
const FORM = `${PAGE} [data-u2="func-form"]`;
const RUN = `${FORM} [data-u2="ff-run"]`;
const ICON = `${FORM} [data-u2="ff-history-icon"]`;

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'func-history'));
  await openApp(page, APP.package, 'u2DemoApp');
  await page.waitForSelector(FI, {timeout: 30000});
}

async function checkPickBuildsForm(page) {
  await page.click(FI);
  await page.waitForSelector('.u2-function-input-popup [data-u2="functions-browser"]', {timeout: 5000});
  await page.fill('.u2-function-input-popup [data-u2="fb-search"] input', 'Sin');
  await page.waitForTimeout(400);
  const picked = await page.evaluate(() => {
    // the fixture qualifies core functions (`core:Sin`); match either spelling
    const row = [...document.querySelectorAll('.u2-function-input-popup [data-u2-func]')]
      .find((r) => /(^|:)Sin$/.test(r.getAttribute('data-u2-func')) && r.offsetParent != null);
    row?.click();
    return row != null;
  });
  await page.waitForSelector(FORM, {timeout: 8000});
  const shape = await page.evaluate(({form, run, icon}) => ({
    editorText: document.querySelector('.u2-function-input-name')?.textContent ?? '',
    inputs: document.querySelectorAll(`${form} .u2-input-root`).length,
    run: document.querySelector(run) != null,
    runInButtons: document.querySelector(`${form} .u2-form-buttons [data-u2="ff-run"]`) != null,
    icon: document.querySelector(icon) != null,
  }), {form: FORM, run: RUN, icon: ICON});
  await shot(page, 'func-call-history-1-form');
  ok('func-call-history/1a picking a function through the FunctionInput builds the form with Run and the history icon',
    picked && shape.inputs > 0 && shape.run && shape.runInButtons && shape.icon,
    JSON.stringify(shape));
}

async function checkRunGate(page) {
  const before = await page.evaluate((run) => {
    const b = document.querySelector(run);
    return {disabled: b?.disabled ?? null, title: b?.title ?? ''};
  }, RUN);
  await page.evaluate((form) => {
    const input = document.querySelector(`${form} .u2-input-root input`);
    input.value = '0.5';
    input.dispatchEvent(new Event('input', {bubbles: true}));
  }, FORM);
  await page.waitForTimeout(300);
  const after = await page.evaluate((run) =>
    document.querySelector(run)?.disabled ?? null, RUN);
  ok('func-call-history/2a Run is gated on the required field and arms when it is filled',
    before.disabled === true && before.title.startsWith('Missing') && after === false,
    `before=${JSON.stringify(before)} after=${after}`);
}

async function checkHistoryPopup(page) {
  await page.click(ICON);
  await page.waitForSelector('.u2-func-form-history-popup [data-u2="func-call-history-browser"]',
    {timeout: 5000});
  await shot(page, 'func-call-history-2-popup');
  await page.keyboard.press('Escape');
  await page.waitForFunction(() =>
    document.querySelector('.u2-func-form-history-popup') == null, undefined, {timeout: 5000});
  const value = await page.evaluate((form) =>
    document.querySelector(`${form} .u2-input-root input`)?.value ?? null, FORM);
  ok('func-call-history/3a the history popup opens from the icon and Esc closes it without a write',
    value === '0.5', `value=${value}`);
}

export const checks = [
  {id: 'func-call-history/1 the FunctionInput pick builds the form with Run and history', run: checkPickBuildsForm},
  {id: 'func-call-history/2 the Run gate follows the required field', run: checkRunGate},
  {id: 'func-call-history/3 the history popup opens and closes without a write', run: checkHistoryPopup},
];
