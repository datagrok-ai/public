import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';


category('FilesView', () => {
  test('Creation', async () => {
    const v = DG.FilesView.create();
    expect(v instanceof DG.FilesView, true);
  });

  test('ShowTree', async () => {
    const v = DG.FilesView.create();
    expect(v.showTree, true);
    v.showTree = false;
    expect($(v.root).find('.d4-tree-view-root').length, 0);
  });

  test('ShowPreview', async () => {
    const v = DG.FilesView.create();
    expect(v.showPreview, true);
    v.showPreview = false;
  });

  test('ShowSearch', async () => {
    const v = DG.FilesView.create();
    grok.shell.addView(v);
    expect(v.showSearch, true);
    v.showSearch = false;
    expect($(v.root).find('.grok-card-view')[0]?.style?.display, 'none');
    v.close();
  });

  test('ShowTreeOnly', async () => {
    const v = DG.FilesView.create();
    expect(v.showTreeOnly, false);
    v.showTreeOnly = true;
    expect(v.showTree, true);
    expect(v.showPreview, false);
    expect(v.showSearch, false);
  });
}, { owner: 'aparamonov@datagrok.ai' });
