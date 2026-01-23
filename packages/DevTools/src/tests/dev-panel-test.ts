import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, before, test, after, delay} from '@datagrok-libraries/test/src/test';
import {EntityType} from '../constants';


category('Dev panel', () => {
  const subs = [];
  const delayDuration = 1000;
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
  });

  test('Dev panel opens', () => new Promise(async (resolve, reject) => {
    if (!entities.length)
      reject('Failed to find entities for the test');
    let stop: boolean = false;
    setTimeout(() => {
      stop = true;
      reject('Timeout exceeded');
    }, delayDuration * entities.length + 100);

    function check() {
      if (devPane != null) {
        stop = false;
        resolve('OK');
      }
      if (!stop)
        setTimeout(check, 50);
    }
    check();
    let result = true;
    let object = null;

    subs.push(grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      currentEntity = acc.context;
      devPane = acc.getPane('Dev');
      result = result && devPane != null && (<DG.Entity>currentEntity).id === object.id;
    }));

    for (const ent of entities) {
      grok.shell.o = ent;
      object = ent;
      await delay(delayDuration);
    }

    if (result)
      resolve('OK');
  }));

  test('Dev panel content', () => new Promise(async (resolve, reject) => {
    if (!entities.length)
      reject('Failed to find entities for the test');

    subs.push(grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      devPane = acc.getPane('Dev');
      devPane.expanded = true;
    }));

    grok.shell.o = entities[0];
    await delay(1500);
    if (!devPane)
      reject('Dev panel not constructed');

    const devPaneContainer = devPane.root.querySelector('.dt-dev-pane-container');
    const snippetSection = devPane.root.querySelector('.dt-snippet-section');
    const textArea = devPane.root.querySelector('.dt-textarea-box');

    if (!devPaneContainer)
      reject('Missing dev pane container');
    else if (!snippetSection)
      reject('Missing snippet section');
    else if (!textArea)
      reject('Missing text area');
    resolve('OK');
  }));

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    grok.shell.windows.showProperties = showPropertyPanel;
  });
});
