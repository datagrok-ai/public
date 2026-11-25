/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
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
import {clickOn, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot inner viewer category order", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-category-order.feature", import.meta.url);
  test("Trellis plot inner viewer category order", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names",""],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names",""],["Viewer Type","Scatter plot"]]);
    await session.step(19, "And user sets \"xColumnName\" inner property of trellis plot viewer to \"RACE\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "RACE"));
    await session.step(20, "And user sets \"yColumnName\" inner property of trellis plot viewer to \"AGE\"", () => setInnerProperty(page, "yColumnName", el("trellis plot viewer"), "AGE"));
    await session.step(21, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(22, "Then the \"cells\" reading of trellis plot viewer should be 2", () => readingIs(page, "cells", el("trellis plot viewer"), 2));
    await run.scenario("X Category Order of the inner scatter plots reorders every cell and stays local", async () => {
      await session.step(25, "When user remembers the \"cell signature F\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F", el("trellis plot viewer")));
      await session.step(26, "And user remembers the \"cell signature M\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M", el("trellis plot viewer")));
      await session.step(27, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(28, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(29, "And user clicks on \"Scatter plot\" tab in context panel", () => clickOn(page, el("\"Scatter plot\" tab in context panel")));
      await session.step(30, "And \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
      await session.step(31, "Then \"X Category Order\" property in context panel should be visible", () => shouldBe(page, el("\"X Category Order\" property in context panel"), "visible"));
      await session.step(32, "When user clicks on button in \"X Category Order\" property in context panel", () => clickOn(page, el("button in \"X Category Order\" property in context panel")));
      await session.step(33, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(34, "And \"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog should be absent", () => shouldBe(page, el("\"APPLY TO ALL VIEWERS\" button in \"Order categories: RACE\" dialog"), "absent"));
      await session.step(35, "When user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(36, "Then the \"cell signature F\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F", el("trellis plot viewer")));
      await session.step(37, "And the \"cell signature M\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M", el("trellis plot viewer")));
      await session.step(38, "And the \"categories of RACE\" reading of filter panel should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "categories of RACE", el("filter panel"), "Asian, Black, Caucasian, Other"));
      await session.step(39, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(40, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(41, "Then the \"cell signature F\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F", el("trellis plot viewer")));
      await session.step(42, "And the \"cell signature M\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature M", el("trellis plot viewer")));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
