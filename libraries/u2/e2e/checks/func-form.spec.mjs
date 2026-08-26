/* W1 of the u2 func-call editor: the U2Demo A/B bench — the
   platform's `DG.InputForm.forFuncCall` and u2's `funcForm` over ONE shared FuncCall
   (packages/U2Demo/src/func-convergence.ts). Boot lands on the tab via the demo app's own
   `u2demo.tab` memory; the page's status lines (`sync =`, `fields =`, `inputs =`, `isValid =`)
   are the bench's own evidence of what the call holds and which side moved. The Dart→u2 leg is
   driven through the Dart input (the same `FuncCallParam.value=` write a `setParamValue` makes —
   the bench exposes no handle to the call itself). */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

/** The bench page, told apart from the inputs-convergence page (same classes, overlapping row
 * names) by the gate row only it has. */
const PAGE = '.u2demo-page:has([data-row="batchId"])';

const u2Cell = (name) => `${PAGE} .u2-input-root[data-row="${name}"]`;
const dartCell = (name) => `${PAGE} .ui-input-root[data-row="${name}"]`;
const u2Input = (name) => `${u2Cell(name)} input`;
const dartInput = (name) => `${dartCell(name)} input.ui-input-editor`;

const statuses = (page) => page.evaluate((sel) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []].map((s) => s.textContent), PAGE);

const status = async (page, prefix) => (await statuses(page)).find((t) => t.startsWith(prefix)) ?? '';

const counters = async (page) => {
  const m = (await status(page, 'sync = ')).match(/u2 → Dart (\d+) · Dart → u2 (\d+)/);
  return m === null ? {toDart: -1, toU2: -1} : {toDart: Number(m[1]), toU2: Number(m[2])};
};

const value = (page, sel) => page.evaluate((sel) => document.querySelector(sel)?.value ?? null, sel);

const waitValue = (page, sel, expected) => page.waitForFunction(({sel, expected}) =>
  document.querySelector(sel)?.value === expected, {sel, expected}, {timeout: 8000});

const waitInputsLine = (page, fragment) => page.waitForFunction(({sel, fragment}) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []]
    .some((s) => s.textContent.startsWith('inputs = ') && s.textContent.includes(fragment)),
{sel: PAGE, fragment}, {timeout: 8000});

/** One macrotask past the last observed state — what an echo storm would need to show itself. */
const settle = (page) => page.evaluate(() => new Promise((resolve) => setTimeout(resolve, 120)));

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'func-convergence'));
  await openApp(page, APP.package, 'u2DemoApp');
  // the grids populate async (forFuncCall awaits); the gate row is appended last
  await page.waitForSelector(`.u2demo-ab [data-row="secret"]`, {timeout: 30000});
  await page.waitForSelector(`.u2demo-ab [data-row="batchId"]`, {timeout: 30000});
}

async function checkBoot(page) {
  const fields = await status(page, 'fields = ');
  const sync = await status(page, 'sync = ');
  const shape = await page.evaluate((sel) => {
    const root = document.querySelector(sel);
    return {
      dart: root.querySelectorAll('.u2demo-ab .ui-input-root').length,
      u2: root.querySelectorAll('.u2demo-ab .u2-input-root').length,
      rows: root.querySelectorAll('.u2demo-ab .u2demo-code[data-row]').length,
      cats: [...root.querySelectorAll('.u2demo-ab-cat')].map((c) => c.textContent),
    };
  }, PAGE);
  await shot(page, 'func-form-1-rest');
  ok('func-form/1a/the-tab-opens-with-both-columns-populated-11-fields-each',
    fields === 'fields = Dart 11 · u2 11 · unsupported: (none)' && shape.dart === shape.u2 && shape.rows === 12 &&
    shape.cats.includes('Advanced'),
    `fields="${fields}" dart roots=${shape.dart} u2 roots=${shape.u2} rows=${shape.rows} ` +
    `cats=${JSON.stringify(shape.cats)}`);

  const inputs = await status(page, 'inputs = ');
  ok('func-form/1b/rest-state-counters-zero-defaults-display-only-u2-shows-them-dart-does-not',
    sync === 'sync = u2 → Dart 0 · Dart → u2 0' &&
    await value(page, u2Input('replicates')) === '3' && await value(page, u2Input('doseLevel')) === '250' &&
    await value(page, dartInput('replicates')) === '' && await value(page, dartInput('doseLevel')) === '' &&
    inputs.includes('replicates: null') && inputs.includes('doseLevel: null'),
    `sync="${sync}" u2 replicates="${await value(page, u2Input('replicates'))}" ` +
    `u2 doseLevel="${await value(page, u2Input('doseLevel'))}" dart replicates=` +
    `"${await value(page, dartInput('replicates'))}" inputs line: ${inputs.slice(0, 200)}`);
}

async function checkU2ToDart(page) {
  const before = await counters(page);
  await page.locator(u2Input('secret')).fill('wo7-u2');
  await waitValue(page, dartInput('secret'), 'wo7-u2');
  await waitInputsLine(page, "secret: 'wo7-u2'");
  const after = await counters(page);
  await settle(page);
  const later = await counters(page);
  ok('func-form/2a/a-u2-edit-reaches-the-dart-input-and-the-call-one-tick-no-echo',
    after.toDart === before.toDart + 1 && after.toU2 === before.toU2 &&
    later.toDart === after.toDart && later.toU2 === after.toU2,
    `dart shows "${await value(page, dartInput('secret'))}"; counters ${JSON.stringify(before)} → ` +
    `${JSON.stringify(after)} → settled ${JSON.stringify(later)}`);
}

async function checkDartToU2(page) {
  const before = await counters(page);
  await page.locator(dartInput('secret')).fill('wo7-dart');
  await waitValue(page, u2Input('secret'), 'wo7-dart');
  await waitInputsLine(page, "secret: 'wo7-dart'");
  const after = await counters(page);
  await settle(page);
  const later = await counters(page);
  ok('func-form/3a/a-dart-edit-reaches-the-u2-input-and-the-call-u2-counter-untouched',
    after.toDart === before.toDart && after.toU2 === before.toU2 + 1 &&
    later.toDart === after.toDart && later.toU2 === after.toU2,
    `u2 shows "${await value(page, u2Input('secret'))}"; counters ${JSON.stringify(before)} → ` +
    `${JSON.stringify(after)} → settled ${JSON.stringify(later)}`);
}

const dartChoiceText = (page) => page.evaluate((sel) => {
  const select = document.querySelector(sel);
  return select === null ? null : (select.selectedOptions[0]?.textContent ?? select.value);
}, `${dartCell('stage')} select`);

async function checkNullableChoice(page) {
  const options = await page.evaluate((sel) =>
    [...document.querySelector(sel)?.options ?? []].map((o) => o.value), `${u2Cell('stage')} select`);
  await page.locator(`${u2Cell('stage')} select`).selectOption('Phase I');
  await waitInputsLine(page, "stage: 'Phase I'");
  const picked = await dartChoiceText(page);
  await page.locator(`${u2Cell('stage')} select`).selectOption('');
  await waitInputsLine(page, 'stage: null');
  const cleared = await dartChoiceText(page);
  ok('func-form/4a/the-nullable-choice-leads-with-an-empty-option-and-a-pick-and-a-clear-round-trip',
    options[0] === '' && options.length === 4 && picked === 'Phase I' && (cleared === '' || cleared === null),
    `u2 options=${JSON.stringify(options)}; picked → dart shows "${picked}", cleared → dart shows "${cleared}"`);
}

const invalidMarks = (page) => page.evaluate(({u2Sel, dartSel}) => ({
  u2: document.querySelector(`${u2Sel} .u2-invalid`) !== null,
  error: document.querySelector(`${u2Sel} .u2-input-error`)?.textContent ?? '',
  dart: document.querySelector(`${dartSel} .d4-invalid`) !== null,
}), {u2Sel: u2Cell('replicates'), dartSel: dartCell('replicates')});

async function checkValidateNotClamp(page) {
  await page.locator(u2Input('replicates')).fill('99');
  await page.waitForFunction(({u2Sel, dartSel}) =>
    document.querySelector(`${u2Sel} .u2-invalid`) !== null && document.querySelector(`${dartSel} .d4-invalid`) !== null,
  {u2Sel: u2Cell('replicates'), dartSel: dartCell('replicates')}, {timeout: 8000});
  const flagged = await invalidMarks(page);
  const kept = await value(page, u2Input('replicates'));
  await shot(page, 'func-form-5-invalid');
  ok('func-form/5a/99-stays-99-and-both-sides-flag-it',
    kept === '99' && flagged.u2 && /at most 10/.test(flagged.error) && flagged.dart,
    `u2 field="${kept}" marks=${JSON.stringify(flagged)}`);

  await page.locator(u2Input('replicates')).fill('5');
  await page.waitForFunction(({u2Sel, dartSel}) =>
    document.querySelector(`${u2Sel} .u2-invalid`) === null && document.querySelector(`${dartSel} .d4-invalid`) === null,
  {u2Sel: u2Cell('replicates'), dartSel: dartCell('replicates')}, {timeout: 8000});
  const cleared = await invalidMarks(page);
  ok('func-form/5b/a-valid-value-clears-both-sides',
    !cleared.u2 && cleared.error === '' && !cleared.dart && await value(page, u2Input('replicates')) === '5',
    JSON.stringify(cleared));
}

const waitGate = (page, expected) => page.waitForFunction(({sel, expected}) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []]
    .some((s) => s.textContent === expected), {sel: PAGE, expected}, {timeout: 8000});

async function checkIsValidGate(page) {
  const rest = await status(page, 'isValid = ');
  await page.locator(`${PAGE} .u2-input-root[data-row="batchId"] input`).fill('B-42');
  await waitGate(page, 'isValid = Dart true · u2 true');
  const filled = await status(page, 'isValid = ');
  await page.locator(`${PAGE} .u2-input-root[data-row="batchId"] input`).fill('');
  await waitGate(page, 'isValid = Dart false · u2 false');
  const cleared = await status(page, 'isValid = ');
  ok('func-form/6a/filling-and-clearing-batchid-from-u2-flips-isvalid-on-both-forms',
    rest === 'isValid = Dart true · u2 false' && filled === 'isValid = Dart true · u2 true' &&
    cleared === 'isValid = Dart false · u2 false',
    `rest="${rest}" (the documented lazy/eager divergence) → filled="${filled}" → cleared="${cleared}"`);
}

async function checkPostfixAndCaption(page) {
  const row = await page.evaluate(({nameSel, dartSel, u2Sel}) => ({
    name: document.querySelector(nameSel)?.textContent ?? '',
    dart: document.querySelector(dartSel)?.textContent ?? '',
    u2: document.querySelector(u2Sel)?.textContent ?? '',
  }), {nameSel: `${PAGE} .u2demo-code[data-row="doseLevel"]`,
    dartSel: dartCell('doseLevel'), u2Sel: u2Cell('doseLevel')});
  ok('func-form/7a/the-doselevel-row-carries-the-camelcased-caption-and-the-mg-postfix-on-both-sides',
    row.name === 'doseLevel' && row.dart.includes('Dose Level') && row.dart.includes('mg') &&
    row.u2.includes('Dose Level') && row.u2.includes('mg'),
    `name cell="${row.name}" dart cell="${row.dart.slice(0, 80)}" u2 cell="${row.u2.slice(0, 80)}"`);
}

export const checks = [
  {id: 'func-form/1 boot: the FuncCall-convergence tab, 11 fields a side, rest state', run: checkBoot},
  {id: 'func-form/2 u2 → Dart: one edit, one counter tick, stable a macrotask later', run: checkU2ToDart},
  {id: 'func-form/3 Dart → u2: the param write follows back, u2 counter untouched', run: checkDartToU2},
  {id: 'func-form/4 nullable choice: empty first option, pick and clear round-trip', run: checkNullableChoice},
  {id: 'func-form/5 validate, not clamp: out-of-range flags both sides and clears', run: checkValidateNotClamp},
  {id: 'func-form/6 isValid gate: batchId filled and cleared from the u2 side', run: checkIsValidGate},
  {id: 'func-form/7 postfix and caption on the doseLevel row', run: checkPostfixAndCaption},
];
