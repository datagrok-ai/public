import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow, withAttachedView} from '../helpers';

// CardView search/filter surface. View must be attached before DOM-backed getters/refresh() work.
// The data source comes from the `dataSource` create option (objectType is only render metadata) —
// without it refresh() throws "DataSource not specified" as an unhandled async rejection.
category('AI: App: CardView Search', () => {
  const withCardView = (body: (cv: DG.CardView) => void): Promise<void> =>
    // awaitRender=false: search/filter state is model-backed; dataSource lets refresh() run cleanly.
    withAttachedView<DG.CardView>(() => {
      const cv = DG.CardView.create({dataSource: grok.dapi.scripts});
      cv.objectType = 'Script';
      return cv;
    }, (cv) => {
      expect(cv instanceof DG.CardView, true);
      body(cv);
    }, false);

  test('searchValue round-trips and survives refresh', async () => {
    await withCardView((cv) => {
      cv.searchValue = 'age > 30';
      expect(cv.searchValue, 'age > 30');
      // refresh() re-runs the search; the value must persist, not be reset.
      expectNoThrow(() => cv.refresh());
      expect(cv.searchValue, 'age > 30');
      cv.searchValue = '';
      expect(cv.searchValue, '');
    });
  });

  test('permanentFilter round-trips', async () => {
    await withCardView((cv) => {
      cv.permanentFilter = '#favorites';
      expect(cv.permanentFilter, '#favorites');
      expectNoThrow(() => cv.refresh());
      expect(cv.permanentFilter, '#favorites');
    });
  });

  test('searchFields round-trips', async () => {
    await withCardView((cv) => {
      cv.searchFields = ['name', 'description'];
      expect(Array.isArray(cv.searchFields), true);
      expect(cv.searchFields.length, 2);
      expect(cv.searchFields.indexOf('name') >= 0, true);
      expect(cv.searchFields.indexOf('description') >= 0, true);
    });
  });

  test('hierarchy + showTree drive grouping state', async () => {
    await withCardView((cv) => {
      cv.hierarchy = ['author'];
      expect(Array.isArray(cv.hierarchy), true);
      expect(cv.hierarchy.indexOf('author') >= 0, true);
      cv.showTree = true;
      expect(cv.showTree, true);
      cv.showTree = false;
      expect(cv.showTree, false);
    });
  });

  test('categoryFilters + filters maps round-trip', async () => {
    await withCardView((cv) => {
      cv.categoryFilters = {tag: 'tags'};
      expect(cv.categoryFilters != null, true);
      expect(cv.categoryFilters['tag'], 'tags');
      cv.filters = {name: 'name'};
      expect(cv.filters != null, true);
      expect(cv.filters['name'], 'name');
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
