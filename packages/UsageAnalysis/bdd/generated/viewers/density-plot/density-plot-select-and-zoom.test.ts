/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/density-plot/density-plot-select-and-zoom.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.density-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, clearSelection, noneSelected, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragZoomOverArea, hasArea, hasNoArea, hoverArea, narrowerRange, noErrors, pickFromContextMenu, readingAtLeast, readingHigher, readingIs, readingLower, readingNotAsRemembered, readingReads, readingsEqual, rememberRange, rememberReading, rememberedRange, repainted, setProperties, setProperty, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Density plot selection, zoom and axis bounds", () => {
  const session = feature(test, "features/viewers/density-plot/density-plot-select-and-zoom.feature", import.meta.url);
  test("Density plot selection, zoom and axis bounds", {tag: ["@journey", "@viewers", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"]]));
    await session.step(20, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
    await session.step(21, "And the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
    await session.step(22, "And density plot viewer should have a \"densest bin\" area", () => hasArea(page, el("density plot viewer"), "densest bin"));
    await run.scenario("Clicking the densest bin selects exactly the rows counted in it", async () => {
      await session.step(25, "Given user clears the row selection", () => clearSelection(page));
      await session.step(26, "Then the \"rows in densest bin\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "rows in densest bin", el("density plot viewer"), 1));
      await session.step(27, "When user clicks on the \"densest bin\" area of density plot viewer", () => clickArea(page, "densest bin", el("density plot viewer")));
      await session.step(28, "Then the \"rows selected\" and \"rows in densest bin\" readings of density plot viewer should be the same", () => readingsEqual(page, "rows selected", "rows in densest bin", el("density plot viewer")));
      await session.step(29, "When user clears the row selection", () => clearSelection(page));
      await session.step(30, "Then no rows should be selected", () => noneSelected(page));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The densest bin's tooltip reports how many rows fell into it", async () => {
      await session.step(34, "When user hovers over the \"densest bin\" area of density plot viewer", () => hoverArea(page, "densest bin", el("density plot viewer")));
      await session.step(35, "Then tooltip should contain text \"rows\"", () => shouldContainText(page, el("tooltip"), "rows"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An Alt-drag zooms into the box and Reset View restores the viewport", async () => {
      await session.step(39, "When user picks \"Reset View\" from the context menu of density plot viewer", () => pickFromContextMenu(page, "Reset View", el("density plot viewer")));
      await session.step(40, "And user remembers the value range of density plot viewer", () => rememberRange(page, el("density plot viewer")));
      await session.step(41, "And user drags a zoom box over the \"view\" area of density plot viewer", () => dragZoomOverArea(page, "view", el("density plot viewer")));
      await session.step(42, "Then density plot viewer should show a narrower value range than before", () => narrowerRange(page, el("density plot viewer")));
      await session.step(43, "And the \"x axis span\" reading of density plot viewer should be lower than before", () => readingLower(page, "x axis span", el("density plot viewer")));
      await session.step(44, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(45, "When user picks \"Reset View\" from the context menu of density plot viewer", () => pickFromContextMenu(page, "Reset View", el("density plot viewer")));
      await session.step(46, "Then density plot viewer should show the remembered value range", () => rememberedRange(page, el("density plot viewer")));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wheel zooms and Reset View undoes it", async () => {
      await session.step(50, "When user picks \"Reset View\" from the context menu of density plot viewer", () => pickFromContextMenu(page, "Reset View", el("density plot viewer")));
      await session.step(51, "And user remembers the value range of density plot viewer", () => rememberRange(page, el("density plot viewer")));
      await session.step(52, "And user scrolls the mouse wheel up over the \"view\" area of density plot viewer", () => wheelOverArea(page, "up", "view", el("density plot viewer")));
      await session.step(53, "Then density plot viewer should show a narrower value range than before", () => narrowerRange(page, el("density plot viewer")));
      await session.step(54, "When user picks \"Reset View\" from the context menu of density plot viewer", () => pickFromContextMenu(page, "Reset View", el("density plot viewer")));
      await session.step(55, "Then density plot viewer should show the remembered value range", () => rememberedRange(page, el("density plot viewer")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The four bound properties pin the viewport, and clearing them releases it", async () => {
      await session.step(59, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["xMin","30"],["xMax","60"]]));
      await session.step(62, "Then the \"x axis min\" reading of density plot viewer should be 30", () => readingIs(page, "x axis min", el("density plot viewer"), 30));
      await session.step(63, "And the \"x axis max\" reading of density plot viewer should be 60", () => readingIs(page, "x axis max", el("density plot viewer"), 60));
      await session.step(64, "And the \"x axis span\" reading of density plot viewer should be 30", () => readingIs(page, "x axis span", el("density plot viewer"), 30));
      await session.step(65, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["yMin","60"],["yMax","100"]]));
      await session.step(68, "Then the \"y axis min\" reading of density plot viewer should be 60", () => readingIs(page, "y axis min", el("density plot viewer"), 60));
      await session.step(69, "And the \"y axis max\" reading of density plot viewer should be 100", () => readingIs(page, "y axis max", el("density plot viewer"), 100));
      await session.step(70, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(71, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["xMin",""],["xMax",""],["yMin",""],["yMax",""]]));
      await session.step(76, "Then the \"x axis span\" reading of density plot viewer should be higher than before", () => readingHigher(page, "x axis span", el("density plot viewer")));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bin To Range re-bins over the viewport rather than the column range", async () => {
      await session.step(80, "When user drags a zoom box over the \"view\" area of density plot viewer", () => dragZoomOverArea(page, "view", el("density plot viewer")));
      await session.step(81, "And user remembers the \"bins drawn\" reading of density plot viewer", () => rememberReading(page, "bins drawn", el("density plot viewer")));
      await session.step(82, "And user sets \"binToRange\" property of density plot viewer to \"true\"", () => setProperty(page, "binToRange", el("density plot viewer"), "true"));
      await session.step(83, "Then the \"bins drawn\" reading of density plot viewer should not be as remembered", () => readingNotAsRemembered(page, "bins drawn", el("density plot viewer")));
      await session.step(84, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(85, "And the \"bin to range\" reading of density plot viewer should be \"true\"", () => readingReads(page, "bin to range", el("density plot viewer"), "true"));
      await session.step(86, "When user sets \"binToRange\" property of density plot viewer to \"false\"", () => setProperty(page, "binToRange", el("density plot viewer"), "false"));
      await session.step(87, "Then the \"bin to range\" reading of density plot viewer should be \"false\"", () => readingReads(page, "bin to range", el("density plot viewer"), "false"));
      await session.step(88, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis over a column with non-positive values shows the axis warning", async () => {
      await session.step(92, "Given user adds a calculated column \"BELOW_ZERO\" with formula \"${WEIGHT} - 200\"", () => addCalculated(page, "BELOW_ZERO", "${WEIGHT} - 200"));
      await session.step(93, "Then the table should have a column \"BELOW_ZERO\"", () => hasColumn(page, "BELOW_ZERO"));
      await session.step(94, "When user sets \"yColumnName\" property of density plot viewer to \"BELOW_ZERO\"", () => setProperty(page, "yColumnName", el("density plot viewer"), "BELOW_ZERO"));
      await session.step(95, "Then density plot viewer should not have a \"y warning\" area", () => hasNoArea(page, el("density plot viewer"), "y warning"));
      await session.step(96, "When user sets \"yAxisType\" property of density plot viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("density plot viewer"), "logarithmic"));
      await session.step(97, "Then density plot viewer should have a \"y warning\" area", () => hasArea(page, el("density plot viewer"), "y warning"));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
      await session.step(99, "When user sets \"yAxisType\" property of density plot viewer to \"linear\"", () => setProperty(page, "yAxisType", el("density plot viewer"), "linear"));
      await session.step(100, "Then density plot viewer should not have a \"y warning\" area", () => hasNoArea(page, el("density plot viewer"), "y warning"));
      await session.step(101, "When user sets \"yColumnName\" property of density plot viewer to \"WEIGHT\"", () => setProperty(page, "yColumnName", el("density plot viewer"), "WEIGHT"));
      await session.step(102, "And user removes \"BELOW_ZERO\" column", () => removeColumn(page, "BELOW_ZERO"));
      await session.step(103, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting the Y axis redraws without moving the value range", async () => {
      await session.step(106, "When user remembers the value range of density plot viewer", () => rememberRange(page, el("density plot viewer")));
      await session.step(107, "And user sets \"invertYAxis\" property of density plot viewer to \"true\"", () => setProperty(page, "invertYAxis", el("density plot viewer"), "true"));
      await session.step(108, "Then density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(109, "And density plot viewer should show the remembered value range", () => rememberedRange(page, el("density plot viewer")));
      await session.step(110, "When user sets \"invertYAxis\" property of density plot viewer to \"false\"", () => setProperty(page, "invertYAxis", el("density plot viewer"), "false"));
      await session.step(111, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
