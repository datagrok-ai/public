import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { category, before, test, after, expect, delay, assure } from '@datagrok-libraries/utils/src/test';
import { EntityType } from '../constants';


category('Dev panel', () => {
  const subs = [];
  const delayDuration = 500;
  const entities: EntityType[] = [];
  const showPropertyPanel = grok.shell.windows.showProperties;
  let currentEntity: EntityType;
  let devPane: DG.AccordionPane;

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
    grok.shell.windows.showProperties = true;
    subs.push(grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      currentEntity = acc.context;
      devPane = acc.getPane('Dev');
    }));
  });

  test('Dev panel opens', async () => {
    if (!entities.length)
      throw 'Failed to find entities for the test';
    for (const ent of entities) {
      grok.shell.o = ent;
      await delay(delayDuration);
      assure.notNull(devPane);
      expect((<DG.Entity>currentEntity).id, (<DG.Entity>ent).id);
      expect(currentEntity.name, ent.name);
    }
  });

  test('Dev panel content', async () => {
    if (!entities.length)
      throw 'Failed to find entities for the test';
    const ent = entities[0];
    grok.shell.o = ent;
    devPane.expanded = true;
    await delay(delayDuration);
    const devPaneContainer = devPane.root.querySelector('.dt-dev-pane-container');
    const snippetSection = devPane.root.querySelector('.dt-snippet-section');
    const textArea = devPane.root.querySelector('.dt-textarea-box');
    [devPaneContainer, snippetSection, textArea].forEach((el) => assure.notNull(el));
  });

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    grok.shell.windows.showProperties = showPropertyPanel;
  });
});
