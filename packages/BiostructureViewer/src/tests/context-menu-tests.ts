import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {category, test, after, before, expect, delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {BiostructureViewerWindowType} from '../utils/context-menu';

import {awaitGrid} from './utils';

import {_package} from '../package-test';

declare const window: BiostructureViewerWindowType;

category('ContextMenu', () => {
  before(async () => {
    if (!window.$biostructureViewer) window.$biostructureViewer = {};
    window.$biostructureViewer.contextMenuError = null;
  });

  after(async () => {
    window.$biostructureViewer!.contextMenuError = null;
  });

  test('RowHeader', async () => {
    const df = await _package.files.readCsv('pdb_id.csv');
    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid, 1000);

    let menu: DG.Menu | null = null;
    view.grid.onContextMenu.subscribe((eventMenu: DG.Menu) => {
      menu = eventMenu;
      return true;
    });
    const gridOverlay = $(view.grid.root).find('canvas').get()[2];
    const cellId1 = view.grid.cell('id', 1);
    const grb = view.grid.root.getBoundingClientRect();
    const cb = cellId1.bounds;
    const cmPE = new PointerEvent('contextmenu', {
      cancelable: true, bubbles: true, view: window, button: 2,
      clientX: grb.left + cb.left - 3, clientY: grb.top + cb.top + 3
    });
    try {
      gridOverlay.dispatchEvent(cmPE);
      // BiostructureViewer:addContextMenu is calling async as package function
      await awaitGrid(view.grid, 1000);
      const err = window.$biostructureViewer?.contextMenuError;
      if (err) {
        // Rethrow error stored on context menu
        throw err;
      }
    } finally {
      if (menu !== null) {
        // Cleanup or hide menu
        // @ts-ignore
        menu.root.style.display = 'none';
      }
    }
  });
});

