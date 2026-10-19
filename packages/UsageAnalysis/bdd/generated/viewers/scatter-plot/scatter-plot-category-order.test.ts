/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {filterToAnyOf, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, hasNoArea, menuDoesNotList, noErrors, pickFromAreaContextMenu, readingIs, rightClickArea, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, areaLiesBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot category order", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-category-order.feature", import.meta.url);
  test("Scatter plot category order", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(21, "Then the \"marker of row 10\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 10", el("scatter plot viewer"), "left", "marker of row 1"));
    await session.step(22, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
    await run.scenario("Row count on the X axis in both directions", async () => {
      await session.step(25, "When user right-clicks on the \"y axis\" area of scatter plot viewer", () => rightClickArea(page, "y axis", el("scatter plot viewer")));
      await session.step(26, "Then the open menu should not list \"Order...\"", () => menuDoesNotList(page, "Order..."));
      await session.step(27, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(28, "And user picks \"Order...\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Order...", "left edge of x axis", el("scatter plot viewer")));
      await session.step(29, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(30, "When user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(31, "Then the \"marker of row 10\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 10", el("scatter plot viewer"), "left", "marker of row 2"));
      await session.step(32, "And the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 1"));
      await session.step(33, "When user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(34, "Then the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
      await session.step(35, "And the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 10\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 10"));
      await session.step(36, "When user selects \"Ascending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(37, "And user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(38, "Then the \"marker of row 10\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 10", el("scatter plot viewer"), "left", "marker of row 1"));
      await session.step(39, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
      await session.step(40, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(41, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row count on a categorical Y axis", async () => {
      await session.step(44, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","AGE"],["Y","RACE"]]), [["X","AGE"],["Y","RACE"]]);
      await session.step(47, "Then the \"marker of row 2\" area of scatter plot viewer should lie above the \"marker of row 1\" area", () => areaLies(page, "marker of row 2", el("scatter plot viewer"), "above", "marker of row 1"));
      await session.step(48, "When user right-clicks on the \"left edge of x axis\" area of scatter plot viewer", () => rightClickArea(page, "left edge of x axis", el("scatter plot viewer")));
      await session.step(49, "Then the open menu should not list \"Order...\"", () => menuDoesNotList(page, "Order..."));
      await session.step(50, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(51, "And user picks \"Order...\" from the context menu of the \"y axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Order...", "y axis", el("scatter plot viewer")));
      await session.step(52, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(53, "Then the \"marker of row 1\" area of scatter plot viewer should lie above the \"marker of row 2\" area", () => areaLies(page, "marker of row 1", el("scatter plot viewer"), "above", "marker of row 2"));
      await session.step(54, "And the \"marker of row 2\" area of scatter plot viewer should lie above the \"marker of row 10\" area", () => areaLies(page, "marker of row 2", el("scatter plot viewer"), "above", "marker of row 10"));
      await session.step(55, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(56, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(57, "Then the \"marker of row 2\" area of scatter plot viewer should lie above the \"marker of row 1\" area", () => areaLies(page, "marker of row 2", el("scatter plot viewer"), "above", "marker of row 1"));
      await session.step(58, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
      await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The order holds while Pack and Zoom by Filter packs the categories left by the filter", async () => {
      await session.step(64, "When user sets \"Zoom And Filter\" property of scatter plot viewer to \"pack and zoom by filter\"", () => setProperty(page, "Zoom And Filter", el("scatter plot viewer"), "pack and zoom by filter"));
      await session.step(65, "And user filters rows where \"RACE\" is one of \"Black, Caucasian, Other\"", () => filterToAnyOf(page, "RACE", "Black, Caucasian, Other"));
      await session.step(66, "Then scatter plot viewer should not have a \"marker of row 10\" area", () => hasNoArea(page, el("scatter plot viewer"), "marker of row 10"));
      await session.step(67, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
      await session.step(68, "When user picks \"Order...\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Order...", "left edge of x axis", el("scatter plot viewer")));
      await session.step(69, "And user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(70, "Then the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 1"));
      await session.step(71, "And the \"x axis max\" reading of scatter plot viewer should be 2.1", () => readingIs(page, "x axis max", el("scatter plot viewer"), 2.1));
      await session.step(72, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(73, "And user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(74, "Then the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 2\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 2"));
      await session.step(75, "When user resets the filter", () => resetFilter(page));
      await session.step(76, "And user sets \"Zoom And Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom And Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(77, "Then the \"marker of row 10\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 10", el("scatter plot viewer"), "left", "marker of row 1"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
