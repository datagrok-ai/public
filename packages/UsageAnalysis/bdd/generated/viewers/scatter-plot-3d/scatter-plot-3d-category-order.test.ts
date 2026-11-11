/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot-3d/scatter-plot-3d-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot3d]
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
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, menuDoesNotList, menuLists, noErrors, pickFromAreaContextMenu, readingReads, rightClickArea, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("3D scatter plot category order", () => {
  const session = feature(test, "features/viewers/scatter-plot-3d/scatter-plot-3d-category-order.feature", import.meta.url);
  test("3D scatter plot category order", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["X","RACE"],["Z","SEX"]]), [["X","RACE"],["Z","SEX"]]);
    await session.step(17, "Then the \"x labels\" reading of 3d scatter plot viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(18, "And the \"z labels\" reading of 3d scatter plot viewer should be \"F, M\"", () => readingReads(page, "z labels", el("3d scatter plot viewer"), "F, M"));
    await run.scenario("X Axis... orders the X labels by row count and by an aggregate, and Default gives them back", async () => {
      await session.step(21, "When user picks \"Order > X Axis...\" from the context menu of the \"view\" area of 3d scatter plot viewer", () => pickFromAreaContextMenu(page, "Order > X Axis...", "view", el("3d scatter plot viewer")));
      await session.step(22, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(23, "When user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(24, "And user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(25, "Then the \"x labels\" reading of 3d scatter plot viewer should be \"Caucasian, Other, Black, Asian\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Caucasian, Other, Black, Asian"));
      await session.step(26, "And the \"z labels\" reading of 3d scatter plot viewer should be \"F, M\"", () => readingReads(page, "z labels", el("3d scatter plot viewer"), "F, M"));
      await session.step(27, "When user selects \"Aggregated value\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Aggregated value", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(28, "And user selects \"AGE\" in \"Column\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "AGE", el("\"Column\" input in \"Order categories: RACE\" dialog")));
      await session.step(29, "Then the \"x labels\" reading of 3d scatter plot viewer should be \"Other, Black, Caucasian, Asian\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Other, Black, Caucasian, Asian"));
      await session.step(30, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(31, "Then the \"x labels\" reading of 3d scatter plot viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Asian, Black, Caucasian, Other"));
      await session.step(32, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(33, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An order on Z leaves X alone, and a numerical axis offers none", async () => {
      await session.step(36, "When user sets \"zCategoryOrder\" property of 3d scatter plot viewer to '{\"manual\":[\"M\",\"F\"]}'", () => setProperty(page, "zCategoryOrder", el("3d scatter plot viewer"), "{\"manual\":[\"M\",\"F\"]}"));
      await session.step(37, "Then the \"z labels\" reading of 3d scatter plot viewer should be \"M, F\"", () => readingReads(page, "z labels", el("3d scatter plot viewer"), "M, F"));
      await session.step(38, "And the \"x labels\" reading of 3d scatter plot viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "x labels", el("3d scatter plot viewer"), "Asian, Black, Caucasian, Other"));
      await session.step(39, "When user sets \"X\" property of 3d scatter plot viewer to \"AGE\"", () => setProperty(page, "X", el("3d scatter plot viewer"), "AGE"));
      await session.step(40, "And user right-clicks on the \"view\" area of 3d scatter plot viewer", () => rightClickArea(page, "view", el("3d scatter plot viewer")));
      await session.step(41, "Then the open menu should list \"Order > Z Axis...\"", () => menuLists(page, "Order > Z Axis..."));
      await session.step(42, "And the open menu should not list \"Order > X Axis...\"", () => menuDoesNotList(page, "Order > X Axis..."));
      await session.step(43, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(44, "And user sets \"zCategoryOrder\" property of 3d scatter plot viewer to \"\"", () => setProperty(page, "zCategoryOrder", el("3d scatter plot viewer"), ""));
      await session.step(45, "And user sets \"X\" property of 3d scatter plot viewer to \"RACE\"", () => setProperty(page, "X", el("3d scatter plot viewer"), "RACE"));
      await session.step(46, "Then the \"z labels\" reading of 3d scatter plot viewer should be \"F, M\"", () => readingReads(page, "z labels", el("3d scatter plot viewer"), "F, M"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
