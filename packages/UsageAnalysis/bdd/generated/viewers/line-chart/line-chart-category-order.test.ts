/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
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
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, onlyOfSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, closeContextMenu, menuDoesNotList, menuLists, noErrors, pickFromAreaContextMenu, readingIs, readingReads, rightClickArea, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, areaLiesBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart category order", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-category-order.feature", import.meta.url);
  test("Line chart category order", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["X","RACE"],["Y","AGE"]]), [["X","RACE"],["Y","AGE"]]);
    await session.step(19, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
    await session.step(20, "And the \"point of row 2\" area of line chart viewer should lie to the left of the \"point of row 3\" area", () => areaLiesBeside(page, "point of row 2", el("line chart viewer"), "left", "point of row 3"));
    await run.scenario("Aggregated value orders the points by the averaged Y column, in both directions", async () => {
      await session.step(23, "When user picks \"Order...\" from the context menu of the \"left edge of x axis\" area of line chart viewer", () => pickFromAreaContextMenu(page, "Order...", "left edge of x axis", el("line chart viewer")));
      await session.step(24, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(25, "When user selects \"Aggregated value\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Aggregated value", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(26, "Then \"Column\" input in \"Order categories: RACE\" dialog should contain text \"AGE\"", () => shouldContainText(page, el("\"Column\" input in \"Order categories: RACE\" dialog"), "AGE"));
      await session.step(27, "And the \"point of row 1\" area of line chart viewer should lie to the left of the \"point of row 3\" area", () => areaLiesBeside(page, "point of row 1", el("line chart viewer"), "left", "point of row 3"));
      await session.step(28, "And the \"point of row 3\" area of line chart viewer should lie to the left of the \"point of row 2\" area", () => areaLiesBeside(page, "point of row 3", el("line chart viewer"), "left", "point of row 2"));
      await session.step(29, "And the \"point of row 2\" area of line chart viewer should lie to the left of the \"point of row 4\" area", () => areaLiesBeside(page, "point of row 2", el("line chart viewer"), "left", "point of row 4"));
      await session.step(30, "When user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(31, "Then the \"point of row 4\" area of line chart viewer should lie to the left of the \"point of row 2\" area", () => areaLiesBeside(page, "point of row 4", el("line chart viewer"), "left", "point of row 2"));
      await session.step(32, "And the \"point of row 3\" area of line chart viewer should lie to the left of the \"point of row 1\" area", () => areaLiesBeside(page, "point of row 3", el("line chart viewer"), "left", "point of row 1"));
      await session.step(33, "When user selects \"Ascending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(34, "And user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(35, "Then the \"point of row 2\" area of line chart viewer should lie to the left of the \"point of row 3\" area", () => areaLiesBeside(page, "point of row 2", el("line chart viewer"), "left", "point of row 3"));
      await session.step(36, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(37, "Then the \"x categories\" reading of line chart viewer should be 4", () => readingIs(page, "x categories", el("line chart viewer"), 4));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A region over an ordered category selects that category's rows", async () => {
      await session.step(41, "When user sets \"xCategoryOrder\" property of line chart viewer to '{\"column\":\"RACE\",\"aggr\":\"values\",\"order\":\"asc\"}'", () => setProperty(page, "xCategoryOrder", el("line chart viewer"), "{\"column\":\"RACE\",\"aggr\":\"values\",\"order\":\"asc\"}"));
      await session.step(42, "And user sets \"annotationRegions\" property of line chart viewer to '[{\"type\":\"area\",\"x\":\"RACE\",\"y\":\"avg(AGE)\",\"header\":\"Third\",\"area\":[[1.5,39],[2.5,39],[2.5,47.5],[1.5,47.5]]}]'", () => setProperty(page, "annotationRegions", el("line chart viewer"), "[{\"type\":\"area\",\"x\":\"RACE\",\"y\":\"avg(AGE)\",\"header\":\"Third\",\"area\":[[1.5,39],[2.5,39],[2.5,47.5],[1.5,47.5]]}]"));
      await session.step(43, "Then the \"point of row 4\" area of line chart viewer should lie inside the \"region Third\" area", () => areaLies(page, "point of row 4", el("line chart viewer"), "inside", "region Third"));
      await session.step(44, "When user clicks on the \"bottom edge of region Third\" area of line chart viewer", () => clickArea(page, "bottom edge of region Third", el("line chart viewer")));
      await session.step(45, "Then only rows where \"RACE\" is \"Other\" should be selected", () => onlyOfSelected(page, "RACE", "Other"));
      await session.step(46, "When user clears the row selection", () => clearSelection(page));
      await session.step(47, "And user sets \"annotationRegions\" property of line chart viewer to \"[]\"", () => setProperty(page, "annotationRegions", el("line chart viewer"), "[]"));
      await session.step(48, "And user sets \"xCategoryOrder\" property of line chart viewer to \"\"", () => setProperty(page, "xCategoryOrder", el("line chart viewer"), ""));
      await session.step(49, "Then the \"point of row 2\" area of line chart viewer should lie to the left of the \"point of row 3\" area", () => areaLiesBeside(page, "point of row 2", el("line chart viewer"), "left", "point of row 3"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numerical X axis offers no Order...", async () => {
      await session.step(53, "When user sets \"X\" property of line chart viewer to \"AGE\"", () => setProperty(page, "X", el("line chart viewer"), "AGE"));
      await session.step(54, "And user right-clicks on the \"left edge of x axis\" area of line chart viewer", () => rightClickArea(page, "left edge of x axis", el("line chart viewer")));
      await session.step(55, "Then the open menu should list \"Invert X Axis\"", () => menuLists(page, "Invert X Axis"));
      await session.step(56, "And the open menu should not list \"Order...\"", () => menuDoesNotList(page, "Order..."));
      await session.step(57, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(58, "And user sets \"X\" property of line chart viewer to \"RACE\"", () => setProperty(page, "X", el("line chart viewer"), "RACE"));
      await session.step(59, "Then the \"point of row 2\" area of line chart viewer should lie to the left of the \"point of row 3\" area", () => areaLiesBeside(page, "point of row 2", el("line chart viewer"), "left", "point of row 3"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
