/* The Locations view's New container (R-d): the ribbon button, disabled until a node is picked, and
   what it opens — the Containers app narrowed to that node's subtree, named after it, with every
   draft it creates starting in that location and its site. Nothing is written: the draft is
   discarded with the view. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {domains, DomainApp} from '@datagrok-libraries/u2/src/dg/index.js';
import {newContainerAt, stockroomLocations} from '../package';

category('Stockroom: locations', () => {
  const site = async (): Promise<string> => {
    const row = await grok.dapi.domains.table('stockroom.location')
      .first({filter: 'site = "Main campus" and kind = "site"'});
    return String(row!.id);
  };

  test('the New container button is in the ribbon, disabled until a location is picked', async () => {
    const view = await stockroomLocations() as DG.View;
    grok.shell.addView(view);
    view.temp['ignoreCloseAll'] = true;
    try {
      await delay(500);
      const add = document.querySelector<HTMLButtonElement>('[data-u2="new-container-button"] button, ' +
        'button[data-u2="new-container-button"]');
      expect(add !== null, true, 'the button is in the view ribbon');
      expect(add!.disabled, true, 'and disabled while the tree has no node selected');
      expect(add!.title, 'Select a location first', 'which says what to do first');
    } finally {
      view.close();
    }
  });

  test('New container opens the Containers app on a draft in the picked location', async () => {
    const id = await site();
    const containers = await domains.table('stockroom.container');
    const view = newContainerAt(containers, id,
      {site: 'Main campus', caption: 'Main campus'}) as DG.View;
    view.temp['ignoreCloseAll'] = true;
    const app = DomainApp.of(view)!;
    try {
      expect(app.listSource.query.value, `location_id under "${id}"`, 'the subtree preset');
      expect(view.name, 'Containers in Main campus', 'the view says which location, not the uuid');
      // the entity source is assigned synchronously; the pristine draft row only lands when
      // its refresh resolves, so the row is what the wait is for
      for (let i = 0; i < 100 && (app.entitySource.value?.currentRow.value ?? null) === null; i++)
        await delay(50);
      expect(app.page.value, 'entity');
      expect(app.entity.value, DomainApp.NEW);
      const draft = app.entitySource.value!.currentRow.value!;
      expect(draft.location_id, id, 'the draft starts in the picked location');
      expect(draft.site, 'Main campus', 'and in its site, which the location picker is narrowed by');
    } finally {
      app.session.discard();
      view.close();
    }
  });
});
