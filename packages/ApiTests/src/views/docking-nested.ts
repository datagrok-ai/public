import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import wu from 'wu';

category('View: DockingNested', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  after(async () => {
    // tv.close();
    // grok.shell.closeTable(df);
  });

  test('left-and-right', async () => {
    const tlmClass = 'test-left-of-main';
    const trmClass = 'test-right-of-main';

    const leftDiv = ui.div('LEFT OF MAIN', {classes: tlmClass});
    const rightDiv = ui.div('RIGHT OF MAIN', {classes: trmClass});

    tv.dockManager.dock(leftDiv, DG.DOCK_TYPE.LEFT, null, 'L-M', 0.2);
    tv.dockManager.dock(rightDiv, DG.DOCK_TYPE.RIGHT, null, 'R-M', 0.2);

    await awaitCheck(() => {
      const resLeftDiv = tv.root.querySelector(`div.${tlmClass}`);
      if (!resLeftDiv) throw new Error('LEFT OF MAIN div not found');
      return isLeftOf(resLeftDiv, tv.grid.root);
    }, 'isLeftOf failed', 100);

    await awaitCheck(() => {
      const resRightDiv = tv.root.querySelector(`div.${trmClass}`);
      if (!resRightDiv) throw new Error('LEFT OF MAIN div not found');
      return isRightOf(resRightDiv, tv.grid.root);
    }, 'isRightOf failed', 100);
  });

  test('top-and-down', async () => {
    const ttmClass = 'test-top-of-main';
    const tdmClass = 'test-down-of-main';

    const topDiv = ui.div('TOP OF MAIN', {classes: ttmClass});
    const downDiv = ui.div('DOWN OF MAIN', {classes: tdmClass});

    tv.dockManager.dock(topDiv, DG.DOCK_TYPE.TOP, null, 'T-M', 0.2);
    tv.dockManager.dock(downDiv, DG.DOCK_TYPE.DOWN, null, 'D-M', 0.2);

    await awaitCheck(() => {
      const resTopDiv = tv.root.querySelector(`div.${ttmClass}`);
      if (!resTopDiv) throw new Error('TOP OF MAIN div not found');
      return isTopOf(resTopDiv, tv.grid.root);
    }, 'isTopOf failed', 100);

    await awaitCheck(() => {
      const resDownDiv = tv.root.querySelector(`div.${tdmClass}`);
      if (!resDownDiv) throw new Error('DOWN OF MAIN div not found');
      return isDownOf(resDownDiv, tv.grid.root);
    }, 'isDownOf failed', 100);
  });

  test('down-of-right', async () => {
    const trmClass = 'test-right-of-main';
    const tdrmClass = 'test-down-of-right-of-main';

    const rightDiv = ui.div('RIGHT OF MAIN', {classes: trmClass});
    const downOfRightDiv = ui.div('DOWN OF RIGHT OF MAIN', {classes: tdrmClass});

    const rightDn = tv.dockManager.dock(rightDiv, DG.DOCK_TYPE.RIGHT, null, 'R-M', 0.2);
    tv.dockManager.dock(downOfRightDiv, DG.DOCK_TYPE.DOWN, rightDn, 'D-R-M', 0.2);

    await awaitCheck(() => {
      const resRightDiv = tv.root.querySelector('div.test-right-of-main');
      if (!resRightDiv) throw new Error('RIGHT OF MAIN div not found');

      const resDownOfRightDiv = tv.root.querySelector('div.test-down-of-right-of-main');
      if (!resDownOfRightDiv) throw new Error('DOWN OF RIGHT OF MAIN div not found');

      return isRightOf(resRightDiv, tv.grid.root) && isRightOf(resDownOfRightDiv, tv.grid.root) &&
        isDownOf(resDownOfRightDiv, resRightDiv);
    }, 'DOWN OF RIGHT OF MAIN failed');
  });

  test('fill-of-down-of-right', async () => {
    const trmClass = 'test-right-of-main';
    const tdrmClass = 'test-down-of-right-of-main';
    const tfdrmClass = 'test-down-of-right-of-main';

    const rightDiv = ui.div('RIGHT OF MAIN', {classes: trmClass});
    const downOfRightDiv = ui.div('DOWN OF RIGHT OF MAIN', {classes: tdrmClass});
    const fillOfDownOfRight = ui.div('FILL OF DOWN OF RIGHT', {classes: tfdrmClass});

    const rightDn = tv.dockManager.dock(rightDiv, DG.DOCK_TYPE.RIGHT, null, 'R-M', 0.4);
    const downOfRightDn = tv.dockManager.dock(downOfRightDiv, DG.DOCK_TYPE.DOWN, rightDn, 'D-R-M', 0.2);
    tv.dockManager.dock(fillOfDownOfRight, DG.DOCK_TYPE.FILL, downOfRightDn, 'F-D-R-M');

    await awaitCheck(() => {
      const resRightDiv = tv.root.querySelector(`div.${trmClass}`);
      if (!resRightDiv) throw new Error('RIGHT OF MAIN div not found');

      const resDownOfRightDiv = tv.root.querySelector(`div.${tdrmClass}`);
      if (!resDownOfRightDiv) throw new Error('DOWN OF RIGHT OF MAIN div not found');

      const resFillOfDownOfRightDiv = tv.root.querySelector(`div.${tfdrmClass}`);
      if (!resFillOfDownOfRightDiv) throw new Error('FILL OF DOWN OF RIGHT OF MAIN div not found');

      return (
        isRightOf(resRightDiv, tv.grid.root) &&
        isRightOf(resDownOfRightDiv, tv.grid.root) &&
        isRightOf(resFillOfDownOfRightDiv, tv.grid.root)
      ) && (
        isDownOf(resDownOfRightDiv, resRightDiv) &&
        isDownOf(resFillOfDownOfRightDiv, resRightDiv)
      );
    }, 'FILL OF DOWN OF RIGHT OF MAIN failed', 100);
  });

  test('fill-of-down-and-fill-of-right', async () => {
    const trmClass = 'test-right-of-main';
    const tfrmClass = 'test-fill-of-right-of-main';
    const tdrmClass = 'test-down-of-right-of-main';
    const tfdrmClass = 'test-down-of-right-of-main';

    const rightDiv = ui.div('RIGHT OF MAIN', {classes: trmClass});
    const fillOfRightDiv = ui.div('FILL OF RIGHT OF MAIN', {classes: tfrmClass});
    const downOfRightDiv = ui.div('DOWN OF RIGHT OF MAIN', {classes: tdrmClass});
    const fillOfDownOfRight = ui.div('FILL OF DOWN OF RIGHT', {classes: tfdrmClass});

    const rightDn = tv.dockManager.dock(rightDiv, DG.DOCK_TYPE.RIGHT, null, 'R-M', 0.4);
    const fillOfRightDn = tv.dockManager.dock(fillOfRightDiv, DG.DOCK_TYPE.FILL, rightDn, 'F-R-M');
    // TODO: rightDn.focus();
    const downOfRightDn = tv.dockManager.dock(downOfRightDiv, DG.DOCK_TYPE.DOWN, rightDn, 'D-R-M', 0.2);
    tv.dockManager.dock(fillOfDownOfRight, DG.DOCK_TYPE.FILL, downOfRightDn, 'F-D-R-M');

    await awaitCheck(() => {
      const resRightDiv = tv.root.querySelector(`div.${trmClass}`);
      if (!resRightDiv) throw new Error('RIGHT OF MAIN div not found');

      const resFillOfRightDiv = tv.root.querySelector(`div.${tfrmClass}`);
      if (!resFillOfRightDiv) throw new Error('FILL OF RIGHT OF MAIN div not found');

      const resDownOfRightDiv = tv.root.querySelector(`div.${tdrmClass}`);
      if (!resDownOfRightDiv) throw new Error('DOWN OF RIGHT OF MAIN div not found');

      const resFillOfDownOfRightDiv = tv.root.querySelector(`div.${tfdrmClass}`);
      if (!resFillOfDownOfRightDiv) throw new Error('FILL OF DOWN OF RIGHT OF MAIN div not found');

      return (
        isRightOf(resRightDiv, tv.grid.root) &&
        isRightOf(resFillOfRightDiv, tv.grid.root) &&
        isRightOf(resDownOfRightDiv, tv.grid.root) &&
        isRightOf(resFillOfDownOfRightDiv, tv.grid.root)
      ) && (
        isDownOf(resDownOfRightDiv, resRightDiv) &&
        isDownOf(resFillOfDownOfRightDiv, resRightDiv) &&
        isDownOf(resFillOfDownOfRightDiv, resFillOfRightDiv)
      );
    }, 'FILL OF DOWN AND FILL OF RIGHT OF MAIN failed', 100);
  });
});

function isLeftOf(tgt: Element, ref: Element): boolean {
  const tgtRect = tgt.getBoundingClientRect();
  const refRect = ref.getBoundingClientRect();

  return (tgtRect.right < refRect.left) &&
    (refRect.top <= tgtRect.top && tgtRect.top <= refRect.bottom);
}

function isRightOf(tgt: Element, ref: Element): boolean {
  const tgtRect = tgt.getBoundingClientRect();
  const refRect = ref.getBoundingClientRect();

  return (tgtRect.left > refRect.left) &&
    (refRect.top <= tgtRect.top && tgtRect.top <= refRect.bottom);
}

function isTopOf(tgt: Element, ref: Element): boolean {
  const tgtRect = tgt.getBoundingClientRect();
  const refRect = ref.getBoundingClientRect();

  return (tgtRect.bottom < refRect.top) &&
    (refRect.left <= tgtRect.left && tgtRect.right <= refRect.right);
}

function isDownOf(tgt: Element, ref: Element): boolean {
  const tgtRect = tgt.getBoundingClientRect();
  const refRect = ref.getBoundingClientRect();

  return (tgtRect.top > refRect.bottom) &&
    (refRect.left <= tgtRect.left && tgtRect.right <= refRect.right);
}
