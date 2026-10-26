/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart, core#974, public#2562]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hoverArea, loadLayout, noErrors, propertyShouldBe, readingIs, reportsNoError, saveLayout, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart regressions from closed issues", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-regressions.feature", import.meta.url);
  test("A vertical bar chart with every row filtered out logs nothing (core#974)", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Orientation","vertical"]]), [["Split","RACE"],["Orientation","vertical"]]);
    await session.step(22, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(23, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(24, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(25, "And bar chart viewer should show 0 rows", () => showsRows(page, el("bar chart viewer"), 0));
    await session.step(26, "And the \"bars\" reading of bar chart viewer should be 0", () => readingIs(page, "bars", el("bar chart viewer"), 0));
    await session.step(27, "When user hovers over the \"view\" area of bar chart viewer", () => hoverArea(page, "view", el("bar chart viewer")));
    await session.step(28, "Then no errors should have been logged", () => noErrors(page));
    await session.step(29, "When user resets the filter", () => resetFilter(page));
    await session.step(30, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(31, "And bar chart viewer should show 1000 rows", () => showsRows(page, el("bar chart viewer"), 1000));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("A datetime Split with Split Map month comes back from a layout with its twelve bars (public#2562)", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@realizes:public#2562"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(37, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","STARTED"],["Split Map","month"]]), [["Split","STARTED"],["Split Map","month"]]);
    await session.step(40, "Then the \"bars\" reading of bar chart viewer should be 12", () => readingIs(page, "bars", el("bar chart viewer"), 12));
    await session.step(41, "And bar chart viewer should have a \"bar January\" area", () => hasArea(page, el("bar chart viewer"), "bar January"));
    await session.step(42, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(43, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
    await session.step(44, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
    await session.step(45, "When user loads the saved layout", () => loadLayout(page));
    await session.step(46, "Then \"Split Map\" property of bar chart viewer should be \"month\"", () => propertyShouldBe(page, "Split Map", el("bar chart viewer"), "month"));
    await session.step(47, "And the \"bars\" reading of bar chart viewer should be 12", () => readingIs(page, "bars", el("bar chart viewer"), 12));
    await session.step(48, "And bar chart viewer should have a \"bar January\" area", () => hasArea(page, el("bar chart viewer"), "bar January"));
    await session.step(49, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
});
