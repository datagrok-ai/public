import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, test, testEvent} from '@datagrok-libraries/test/src/test';

class StatusViewer extends DG.JsViewer {
  getWidgetStatus(): DG.IWidgetStatus {
    const status = super.getWidgetStatus();
    return {...status, parts: {...status.parts, content: this.root}, values: {...status.values, original: 5}};
  }
}

function checkNamedProviders(widget: DG.Widget, lookupByRoot: boolean = true): void {
  let count = 2;
  let calls = 0;
  const marker = ui.divText('marker');
  widget.addStatusProvider('renderer', () => {
    calls++;
    return {parts: {marker}, hitAreas: {marker: {x: 11, y: 12, width: 13, height: 14}},
      values: {'custom count': count}};
  });
  const first = widget.getWidgetStatus();
  expect(calls, 1);
  expect(first.parts.marker, marker);
  expect(first.hitAreas.marker.x, 11);
  expect(first.hitAreas.marker.width, 13);
  expect(first.values?.['custom count'], 2);

  count = 4;
  const reader = lookupByRoot ? DG.Widget.find(widget.root)! : widget;
  expect(reader != null, true);
  expect(reader.getWidgetStatus().values?.['custom count'], 4);
  expect(calls, 2);
  widget.addStatusProvider('override', () => ({values: {'custom count': 7}}));
  let replacementCalls = 0;
  widget.addStatusProvider('renderer', () => {
    replacementCalls++;
    return {values: {'custom count': 9}};
  });
  expect(widget.getWidgetStatus().values?.['custom count'], 7);
  expect(calls, 2);
  expect(replacementCalls, 1);
  widget.removeStatusProvider('override');
  expect(widget.getWidgetStatus().values?.['custom count'], 9);
  expect(replacementCalls, 2);
  widget.removeStatusProvider('renderer');
  widget.removeStatusProvider('absent');
  const removed = widget.getWidgetStatus();
  expect(removed.values?.['custom count'] === undefined, true);
  expect(removed.parts.marker === undefined, true);
  expect(removed.hitAreas.marker === undefined, true);
  expect(replacementCalls, 2);

  let failureCalls = 0;
  widget.addStatusProvider('failing', () => {
    failureCalls++;
    throw new Error('status provider failure');
  });
  let propagated = false;
  try {
    widget.getWidgetStatus();
  } catch (error) {
    propagated = String(error).includes('status provider failure');
  }
  expect(propagated, true);
  expect(failureCalls, 1);
  widget.removeStatusProvider('failing');
  widget.getWidgetStatus();
  expect(failureCalls, 1);
}

category('Widget: status providers', () => {
  test('plain Widget merges providers and clears them on detach', async () => {
    const widget = new DG.Widget(ui.div());
    try {
      checkNamedProviders(widget);
      widget.addStatusProvider('detach', () => ({values: {detached: true}}));
      widget.detach();
      expect(widget.getWidgetStatus().values?.detached === undefined, true);
    } finally {
      widget.detach();
    }
  });

  test('JsViewer override composes with super and closing its view clears providers', async () => {
    const viewer = new StatusViewer();
    const view = grok.shell.addTableView(grok.data.demo.demog(10));
    try {
      view.addViewer(viewer);
      checkNamedProviders(viewer, false);
      expect(viewer.getWidgetStatus().values?.original, 5);
      viewer.addStatusProvider('detach', () => ({values: {detached: true}}));
      view.close();
      expect(viewer.getWidgetStatus().values?.detached === undefined, true);
    } finally {
      if (Array.from(grok.shell.tableViews).some((v) => v.id === view.id))
        view.close();
    }
  });

  test('native View aliases share providers and closing the view clears them', async () => {
    const view = grok.shell.newView('Native view status');
    const alias = new DG.View(view.dart);
    try {
      checkNamedProviders(view);
      view.addStatusProvider('shared', () => ({values: {shared: true}}));
      expect(alias.getWidgetStatus().values?.shared, true);
      alias.removeStatusProvider('shared');
      expect(view.getWidgetStatus().values?.shared === undefined, true);
      alias.addStatusProvider('detach', () => ({values: {detached: true}}));
      expect(view.getWidgetStatus().values?.detached, true);
      view.close();
      expect(alias.getWidgetStatus().values?.detached === undefined, true);
    } finally {
      if (Array.from(grok.shell.views).some((v) => v.id === view.id))
        view.close();
    }
  });

  test('native ColumnGrid aliases share providers and dispose with their host', async () => {
    const columns = DG.ColumnGrid.columnSelector(grok.data.demo.demog(10));
    const alias = new DG.ColumnGrid(columns.dart);
    const view = grok.shell.newView('Column grid status', [columns.root]);
    try {
      checkNamedProviders(columns);
      columns.addStatusProvider('shared', () => ({values: {shared: true}}));
      expect(alias.getWidgetStatus().values?.shared, true);
      alias.removeStatusProvider('shared');
      expect(columns.getWidgetStatus().values?.shared === undefined, true);
      alias.addStatusProvider('detach', () => ({values: {detached: true}}));
      view.close();
      expect(columns.getWidgetStatus().values?.detached === undefined, true);
    } finally {
      if (Array.from(grok.shell.views).some((v) => v.id === view.id))
        view.close();
    }
  });

  test('DartWidget providers are shared between wrappers', async () => {
    const accordion = ui.accordion();
    const view = grok.shell.newView('Status providers', [accordion.root]);
    try {
      checkNamedProviders(accordion);
      accordion.addStatusProvider('shared', () => ({values: {shared: true}}));
      const wrapper = new DG.DartWidget(accordion.dart);
      expect(wrapper.getWidgetStatus().values?.shared, true);
      wrapper.removeStatusProvider('shared');
      expect(accordion.getWidgetStatus().values?.shared === undefined, true);
      // detaching a wrapper leaves the native widget, and its providers, alive
      accordion.addStatusProvider('detach', () => ({values: {detached: true}}));
      wrapper.detach();
      expect(accordion.getWidgetStatus().values?.detached, true);
      view.close();
      expect(accordion.getWidgetStatus().values?.detached === undefined, true);
    } finally {
      if (Array.from(grok.shell.views).some((v) => v.id === view.id))
        view.close();
    }
  });

  test('Grid keeps native fields and geometry, and clears providers on rebinding and disposal', async () => {
    const view = grok.shell.addTableView(grok.data.demo.demog(10));
    const grid = view.grid;
    try {
      await testEvent(grid.onAfterDrawContent, () => {}, () => grid.invalidate(), 10000,
        'Grid did not draw its native hit areas');
      checkNamedProviders(grid);
      expect(grid.getWidgetStatus().parts.canvas, grid.canvas);
      expect(grid.getWidgetStatus().values?.rows, 10);
      const area = Object.values(grid.getWidgetStatus().hitAreas)[0];
      expect(typeof area.x, 'number');
      expect(typeof area.width, 'number');
      let oldFrameCalls = 0;
      grid.addStatusProvider('old-frame', () => {
        oldFrameCalls++;
        return {values: {'old frame': true}};
      });
      expect(grid.getWidgetStatus().values?.['old frame'], true);
      grid.dataFrame = grok.data.demo.demog(4);
      expect(grid.getWidgetStatus().values?.['old frame'] === undefined, true);
      expect(oldFrameCalls, 1);
      grid.addStatusProvider('detach', () => ({values: {'removed on detach': true}}));
      view.close();
      expect(grid.getWidgetStatus().values?.['removed on detach'] === undefined, true);
    } finally {
      if (Array.from(grok.shell.tableViews).some((v) => v.id === view.id))
        view.close();
    }
  });
});
