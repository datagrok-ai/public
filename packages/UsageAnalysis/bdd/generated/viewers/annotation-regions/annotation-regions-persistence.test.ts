/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/annotation-regions/annotation-regions-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, loadLayout, noErrors, readingIs, saveLayout, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Annotation regions persist with the view", () => {
  const session = feature(test, "features/viewers/annotation-regions/annotation-regions-persistence.feature", import.meta.url);
  test("The regions and their titles survive a layout round trip", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]);
    await session.step(16, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
    await session.step(17, "And scatter plot viewer should have a \"region Outer title\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer title"));
    await session.step(18, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(19, "And user loads the saved layout", () => loadLayout(page));
    await session.step(20, "Then the \"viewer regions\" reading of scatter plot viewer should be 2", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 2));
    await session.step(21, "And the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
    await session.step(22, "And scatter plot viewer should have a \"region Outer\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer"));
    await session.step(23, "And scatter plot viewer should have a \"region Outer title\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer title"));
    await session.step(24, "And scatter plot viewer should have a \"region Inner title\" area", () => hasArea(page, el("scatter plot viewer"), "region Inner title"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
  test("The regions survive a layout saved through the server", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]}]"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]}]"]]);
    await session.step(32, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
    await session.step(33, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(34, "And user loads the saved layout", () => loadLayout(page));
    await session.step(35, "Then the \"viewer regions\" reading of scatter plot viewer should be 1", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 1));
    await session.step(36, "And scatter plot viewer should have a \"region Outer\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer"));
    await session.step(37, "And scatter plot viewer should have a \"region Outer title\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer title"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
});
