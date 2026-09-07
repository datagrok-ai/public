/* Card and StatCard as spec nodes: both tags render from a dg-ui/1 document, the card takes spec
   children into its body, a click in Run mode toggles the selection ring, and a stat card bound
   to a `u2-state` follows the value an input writes into it. */
import {ok, shot} from '../local.mjs';
import {expandTree, openSpec, reopenApp, setField, specErrors, surfaceNode, toMode,
  waitStatus} from '../lib.mjs';

const CARD_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-state', name: 'metric', props: {type: 'string', initial: '1.2M'}}],
  root: {tag: 'u2-panel', name: 'cardRoot', children: [
    {tag: 'u2-card', name: 'compound',
      props: {title: 'Aspirin', subtitle: 'NSAID', icon: 'capsules', selectable: true},
      children: [{tag: 'u2-text-input', name: 'metricInput', props: {label: 'Metric'},
        bind: {value: '$.metric'}}]},
    {tag: 'u2-stat-card', name: 'revenue', props: {label: 'Revenue', delta: 0.12},
      bind: {value: '$.metric'}},
  ]},
}, null, 2);

export async function fixture(page) {
  // the file before this one leaves the app on another view: the designer ribbon has to be back
  // before its Open menu can be walked
  await reopenApp(page);
  await openSpec(page, CARD_SPEC);
  await waitStatus(page, 'cardRoot');
  await expandTree(page);
}

const partText = (page, name, cls) => page.locator(`.u2-designer-surface [data-u2-name="${name}"] ${cls}`)
  .first().textContent();

async function checkRendering(page) {
  const cards = await surfaceNode(page, 'compound').count();
  const stats = await surfaceNode(page, 'revenue').count();
  const errors = await specErrors(page);
  const title = await partText(page, 'compound', '.u2-card-title');
  const subtitle = await partText(page, 'compound', '.u2-card-subtitle');
  const label = await partText(page, 'revenue', '.u2-stat-label');
  await shot(page, 'card-1-rendered');
  ok('card/1a/both-tags-render-placeholder-free',
    cards === 1 && stats === 1 && errors === 0 && title === 'Aspirin' && subtitle === 'NSAID' &&
    label === 'Revenue',
    `cards=${cards} stats=${stats} errors=${errors} title="${title}" ` +
    `subtitle="${subtitle}" label="${label}"`);

  const inBody = await page.locator(
    '.u2-designer-surface [data-u2-name="compound"] .u2-card-body [data-u2-name="metricInput"]').count();
  ok('card/1b/spec-children-land-in-the-card-body', inBody === 1, `inputs in body=${inBody}`);

  const value = await partText(page, 'revenue', '.u2-stat-value');
  const delta = await partText(page, 'revenue', '.u2-stat-delta');
  ok('card/1c/the-stat-card-shows-the-bound-value-and-its-delta',
    value === '1.2M' && delta === '+12%', `value="${value}" delta="${delta}"`);
}

async function checkInteraction(page) {
  await toMode(page, 'Run');
  const card = surfaceNode(page, 'compound');
  const before = await card.evaluate((el) => el.classList.contains('u2-card-selected'));
  await card.locator('.u2-card-title').click();
  await page.waitForTimeout(300);
  const after = await card.evaluate((el) => el.classList.contains('u2-card-selected'));
  await shot(page, 'card-2-selected');
  ok('card/2a/a-click-toggles-the-selection-ring', before === false && after === true,
    `before=${before} after=${after}`);

  await setField(page, 'metricInput', '2.4M');
  const value = await partText(page, 'revenue', '.u2-stat-value');
  await shot(page, 'card-3-bound-value');
  ok('card/2b/the-stat-card-follows-the-state-an-input-writes', value === '2.4M', `value="${value}"`);
  await toMode(page, 'Design');
}

export const checks = [
  {id: 'card/1 card and stat card render from a spec', run: checkRendering},
  {id: 'card/2 click selection and the bound stat value', run: checkInteraction},
];
