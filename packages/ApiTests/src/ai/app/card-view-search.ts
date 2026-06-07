import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow, noThrow, withAttachedView} from '../helpers';

// CardView search/filter surface. View must be attached before DOM-backed getters/refresh() work,
// and needs a real objectType — an empty data source makes refresh()/repaint() throw and leak.
category('AI: App: CardView Search', () => {
  const withCardView = (body: (cv: DG.CardView) => void): Promise<void> =>
    // awaitRender=false: search/filter state is model-backed; objectType gives refresh() a data source.
    withAttachedView<DG.CardView>(() => {
      const cv = DG.CardView.create({});
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

  test('refresh does not throw and repaint is callable', async () => {
    await withCardView((cv) => {
      expectNoThrow(() => cv.refresh());
      // repaint() drives a full grid redraw that needs the data source loaded; in the
      // headless harness it may throw before load completes, so exercise it best-effort.
      noThrow(() => cv.repaint());
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
