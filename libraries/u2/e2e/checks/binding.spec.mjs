/* Universal binding (UB-4/UB-10) — the two live tiers as a user feels them: an h3's text bound
   through the Properties row button follows a Run-mode edit with the element never re-created,
   and an input's label follows per keystroke while the focus stays in the sibling field being
   typed into. */
import {ok, shot} from '../local.mjs';
import {bindThroughPicker, clearBalloons, expandTree, openSpec, pickLeaf, reopenApp, selectRow,
  setField, toMode, waitStatus} from '../lib.mjs';

const BIND_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-state', name: 'who', props: {type: 'string', initial: 'world'}}],
  root: {tag: 'u2-form', name: 'bindForm', children: [
    {tag: 'h3', name: 'title', props: {text: 'Untitled'}},
    {tag: 'u2-text-input', name: 'driver', props: {label: 'Driver'}},
    {tag: 'u2-text-input', name: 'sibling', props: {label: 'Sibling'}},
  ]},
}, null, 2);

export async function fixture(page) {
  // the file before this one leaves the app on another view with error balloons still up: the
  // designer ribbon has to be back, and the balloons off the panel they cover
  await reopenApp(page);
  await clearBalloons(page);
  await openSpec(page, BIND_SPEC);
  await waitStatus(page, 'bindForm');
  await expandTree(page);
}

const titleState = (page) => page.evaluate(() => {
  const el = document.querySelector('.u2-designer-surface [data-u2-name="title"]');
  return {text: el?.textContent ?? null, probe: el?.__probe ?? null};
});

async function checkHtmlLiveText(page) {
  await selectRow(page, 'title');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="text"]').first().click();
  const picked = await pickLeaf(page, 'value', 'who');
  await bindThroughPicker(page, 'driver', 'value', 'who');
  const design = await titleState(page);
  await toMode(page, 'Run');
  // stamped after entering Run: the mode toggle itself rebuilds every node bound to a tray
  // component (DD9) — what must be identity-stable is the Run-mode write
  await page.evaluate(() =>
    document.querySelector('.u2-designer-surface [data-u2-name="title"]').__probe = 1);
  await setField(page, 'driver', 'hello');
  const run = await titleState(page);
  await shot(page, 'binding-1-html-live');
  await toMode(page, 'Design');
  ok('binding/1a/a-bound-h3-follows-the-state-with-its-element-never-re-created',
    picked.found !== undefined && design.text === 'world' && run.text === 'hello' && run.probe === 1,
    `picked=${picked.found} design=${JSON.stringify(design)} run=${JSON.stringify(run)}`);
}

async function checkInputLabelLive(page) {
  await bindThroughPicker(page, 'sibling', 'value', 'who');
  await selectRow(page, 'driver');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="label"]').first().click();
  const picked = await pickLeaf(page, 'value', 'who');
  await toMode(page, 'Run');
  const editor = page.locator('.u2-designer-surface [data-u2-name="sibling"] input').first();
  await editor.fill('');
  await editor.click();
  const steps = [];
  for (const ch of 'Live') {
    await page.keyboard.type(ch);
    await page.waitForTimeout(150);
    steps.push(await page.evaluate(() => {
      const input = document.querySelector('.u2-designer-surface [data-u2-name="sibling"] input');
      return {value: input?.value ?? null, focused: document.activeElement === input,
        label: document.querySelector('.u2-designer-surface [data-u2-name="driver"] .u2-input-label')
          ?.textContent ?? null};
    }));
  }
  await shot(page, 'binding-2-label-live');
  await toMode(page, 'Design');
  ok('binding/2a/a-bound-label-follows-per-keystroke-and-the-focus-stays-in-the-sibling',
    picked.found !== undefined && steps.length === 4 &&
    steps.every((s) => s.label === s.value && s.value !== '' && s.focused),
    `picked=${picked.found} steps=${JSON.stringify(steps)}`);
}

export const checks = [
  {id: 'binding/1 an HTML text bind is live: the canvas follows, element identity stable', run: checkHtmlLiveText},
  {id: 'binding/2 an input label bind is live: per-keystroke follow, focus retained', run: checkInputLabelLive},
];
