import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, before, test, after, expect, waitUntil} from '@datagrok-libraries/test/src/test';
import {EntityType} from '../constants';


category('Dev panel', () => {
  const entities: EntityType[] = [];
  const showPropertyPanel = grok.shell.windows.showContextPanel;

  before(async () => {
    const df = grok.data.demo.demog();
    entities.push(...[
      await grok.dapi.connections.first(),
      await grok.dapi.queries.first(),
      await grok.dapi.users.first(),
      await grok.dapi.groups.first(),
      await grok.dapi.packages.first(),
      await grok.dapi.projects.first(),
      await grok.dapi.scripts.first(),
      await grok.dapi.layouts.first(),
      DG.Func.find()[0],
      df,
      df.col('age'),
      DG.View.create(),
    ].filter((ent) => ent != null));
    grok.shell.windows.showContextPanel = true;
  });

  /** Makes [ent] the current object and resolves with the Dev pane the context panel builds for it. */
  async function devPaneOf(ent: EntityType): Promise<DG.AccordionPane> {
    let pane: DG.AccordionPane | null = null;
    const sub = grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      if ((acc.context as DG.Entity)?.id === (ent as DG.Entity).id)
        pane = acc.getPane('Dev');
    });
    try {
      grok.shell.o = ent;
      await waitUntil(() => pane != null, `no Dev pane for ${ent}`);
    } finally {
      sub.unsubscribe();
    }
    return pane!;
  }

  test('Dev panel opens', async () => {
    expect(entities.length > 0, true, 'Failed to find entities for the test');
    for (const ent of entities)
      await devPaneOf(ent);
  });

  test('Dev panel content', async () => {
    expect(entities.length > 0, true, 'Failed to find entities for the test');
    const pane = await devPaneOf(entities[0]);
    pane.expanded = true;
    await waitUntil(() => pane.root.querySelector('.dt-dev-pane-container') != null, 'Missing dev pane container');
    expect(pane.root.querySelector('.dt-snippet-section') != null, true, 'Missing snippet section');
    expect(pane.root.querySelector('.dt-textarea-box') != null, true, 'Missing text area');
  });

  after(async () => {
    grok.shell.windows.showContextPanel = showPropertyPanel;
  });
});
