/* The front door: what u2 is, what the six demo areas cover, and the three sub-demos worth
   reading first. The default leaf when nothing is remembered and no path is in the URL. */
import {divV, divH, span, h3, link} from '@datagrok-libraries/u2';
import {DEMO_TREE, DemoContext, DemoLeaf} from '../nav';

const START_HERE = ['all-inputs', 'basic-inputs', 'form'];

export function overviewPage(ctx: DemoContext): HTMLElement {
  const starts: DemoLeaf[] = [];
  for (const group of DEMO_TREE) {
    for (const leaf of group.children)
      if (START_HERE.includes(leaf.id))
        starts.push(leaf);
  }
  // the group this page itself lives in is not one of the areas it describes
  const areas = DEMO_TREE.filter((g) => g.children.some((l) => l.id !== 'overview'));

  return divV([
    span('u2 is Datagrok\'s UI library. A component is a plain TypeScript class with a DOM root ' +
      'and signal-valued properties: you build it with new, read and write its state through ' +
      'signals, and skin it with design tokens. Nothing on the pages below needs a viewer or a ' +
      'Dart form — the u2/dg layer is what adds the platform.', 'u2demo-hint'),
    span('Every sub-demo shows the source of the build you are running in the context panel, and ' +
      'the ribbon\'s crosshairs toggle turns any control on the page into a live property grid ' +
      'there. Each one has its own URL, so a link reproduces exactly what you are looking at.',
    'u2demo-hint'),
    h3('Start here'),
    divV(starts.map((leaf) => divH([
      link(`${leaf.group.label} / ${leaf.label}`, () => ctx.navigate(leaf.id)),
      span(leaf.description, 'u2demo-dim'),
    ], 'u2demo-row'))),
    h3('What is where'),
    divV(areas.map((group) => divH([
      span(group.label, 'u2demo-overview-area'),
      span(group.children.map((leaf) => leaf.label).join(' · '), 'u2demo-dim'),
    ], 'u2demo-row'))),
  ], 'u2demo-page');
}
