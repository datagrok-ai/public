import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectBoolGetSet, withAttachedView} from '../helpers';

// DG.FilesView (extends CardView): show*/tree getters read the rendered DOM, so unattached
// getters throw on a null context — withAttachedView gives each test a freshly-attached view.
category('AI: App: FilesView JS API', () => {
  const withFilesView = (body: (fv: DG.FilesView) => void | Promise<void>): Promise<void> =>
    withAttachedView<DG.FilesView>(() => DG.FilesView.create({}), body, false);

  test('create() returns a FilesView instance', async () => {
    await withFilesView((fv) => expect(fv instanceof DG.FilesView, true));
  });

  test('showPreview round-trips true/false', async () => {
    await withFilesView((fv) => expectBoolGetSet(fv, 'showPreview'));
  });

  test('showTree round-trips true/false', async () => {
    await withFilesView((fv) => expectBoolGetSet(fv, 'showTree'));
  });

  test('showSearch round-trips true/false', async () => {
    await withFilesView((fv) => expectBoolGetSet(fv, 'showSearch'));
  });

  test('tree getter returns a non-null TreeViewGroup', async () => {
    await withFilesView((fv) => {
      const tree = fv.tree;
      expect(tree != null, true);
      expect(tree instanceof DG.TreeViewGroup, true);
    });
  });

  test('showTreeOnly=true collapses to tree-only configuration', async () => {
    await withFilesView((fv) => {
      fv.showTreeOnly = true;
      expect(fv.showTreeOnly, true);
      // Computed from the three underlying flags.
      expect(fv.showTree, true);
      expect(fv.showPreview, false);
      expect(fv.showSearch, false);
    });
  });

  test('showTreeOnly=false restores preview and search', async () => {
    await withFilesView((fv) => {
      fv.showTreeOnly = true;
      expect(fv.showTreeOnly, true);
      fv.showTreeOnly = false;
      expect(fv.showTreeOnly, false);
      // Setter forces showTree=true and re-enables preview/search.
      expect(fv.showTree, true);
      expect(fv.showPreview, true);
      expect(fv.showSearch, true);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
