/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/density-plot/density-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.density-plot]
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
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, boundTable, hasArea, hasNoArea, lessInk, noErrors, painted, readingAtLeast, readingDiffers, readingIs, readingLower, readingReads, readingsEqual, repainted, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {descriptionAbove, descriptionBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Density plot chrome, binding and the viewer filter", () => {
  const session = feature(test, "features/viewers/density-plot/density-plot.feature", import.meta.url);
  test("Density plot chrome, binding and the viewer filter", {tag: ["@journey", "@viewers", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","HEIGHT"]]));
    await session.step(21, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(22, "And the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
    await session.step(23, "And the \"x column\" reading of density plot viewer should be \"AGE\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE"));
    await session.step(24, "And the \"y column\" reading of density plot viewer should be \"HEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "HEIGHT"));
    await session.step(25, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
    await run.scenario("A blank in either column keeps the row out of the bins", async () => {
      await session.step(28, "Then the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(29, "When user sets \"yColumnName\" property of density plot viewer to \"WEIGHT\"", () => setProperty(page, "yColumnName", el("density plot viewer"), "WEIGHT"));
      await session.step(30, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(31, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(32, "When user sets \"yColumnName\" property of density plot viewer to \"HEIGHT\"", () => setProperty(page, "yColumnName", el("density plot viewer"), "HEIGHT"));
      await session.step(33, "Then the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Color Scale decides whether the scale is drawn at all", async () => {
      await session.step(37, "Then density plot viewer should have a \"color scale\" area", () => hasArea(page, el("density plot viewer"), "color scale"));
      await session.step(38, "And the \"color scale\" area of density plot viewer should be painted", () => areaPainted(page, "color scale", el("density plot viewer")));
      await session.step(39, "When user sets \"showColorScale\" property of density plot viewer to \"false\"", () => setProperty(page, "showColorScale", el("density plot viewer"), "false"));
      await session.step(40, "Then density plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("density plot viewer"), "color scale"));
      await session.step(41, "And density plot viewer should have less ink than before", () => lessInk(page, el("density plot viewer")));
      await session.step(42, "When user sets \"showColorScale\" property of density plot viewer to \"true\"", () => setProperty(page, "showColorScale", el("density plot viewer"), "true"));
      await session.step(43, "Then density plot viewer should have a \"color scale\" area", () => hasArea(page, el("density plot viewer"), "color scale"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The colour scale bounds are reported whether or not the scale is drawn", async () => {
      await session.step(47, "Then the \"color scale min\" reading of density plot viewer should be 0", () => readingIs(page, "color scale min", el("density plot viewer"), 0));
      await session.step(48, "And the \"color scale max\" and \"max bin count\" readings of density plot viewer should be the same", () => readingsEqual(page, "color scale max", "max bin count", el("density plot viewer")));
      await session.step(49, "When user sets \"showColorScale\" property of density plot viewer to \"false\"", () => setProperty(page, "showColorScale", el("density plot viewer"), "false"));
      await session.step(50, "Then the \"color scale max\" and \"max bin count\" readings of density plot viewer should be the same", () => readingsEqual(page, "color scale max", "max bin count", el("density plot viewer")));
      await session.step(51, "When user sets \"showColorScale\" property of density plot viewer to \"true\"", () => setProperty(page, "showColorScale", el("density plot viewer"), "true"));
      await session.step(52, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show X Axis and Show Y Axis take the axis boxes away", async () => {
      await session.step(55, "Then density plot viewer should have an \"x axis\" area", () => hasArea(page, el("density plot viewer"), "x axis"));
      await session.step(56, "And density plot viewer should have a \"y axis\" area", () => hasArea(page, el("density plot viewer"), "y axis"));
      await session.step(57, "When user sets \"showXAxis\" property of density plot viewer to \"false\"", () => setProperty(page, "showXAxis", el("density plot viewer"), "false"));
      await session.step(58, "Then density plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("density plot viewer"), "x axis"));
      await session.step(59, "And density plot viewer should have a \"y axis\" area", () => hasArea(page, el("density plot viewer"), "y axis"));
      await session.step(60, "When user sets \"showYAxis\" property of density plot viewer to \"false\"", () => setProperty(page, "showYAxis", el("density plot viewer"), "false"));
      await session.step(61, "Then density plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("density plot viewer"), "y axis"));
      await session.step(62, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["showXAxis","true"],["showYAxis","true"]]));
      await session.step(65, "Then density plot viewer should have an \"x axis\" area", () => hasArea(page, el("density plot viewer"), "x axis"));
      await session.step(66, "And density plot viewer should have a \"y axis\" area", () => hasArea(page, el("density plot viewer"), "y axis"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show X Selector and Show Y Selector are read as the layout resolved them", async () => {
      await session.step(70, "Then the \"x selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "x selector shown", el("density plot viewer"), "true"));
      await session.step(71, "And the \"y selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "y selector shown", el("density plot viewer"), "true"));
      await session.step(72, "When user sets \"showXSelector\" property of density plot viewer to \"false\"", () => setProperty(page, "showXSelector", el("density plot viewer"), "false"));
      await session.step(73, "Then the \"x selector shown\" reading of density plot viewer should be \"false\"", () => readingReads(page, "x selector shown", el("density plot viewer"), "false"));
      await session.step(74, "And the \"y selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "y selector shown", el("density plot viewer"), "true"));
      await session.step(75, "When user sets \"showYSelector\" property of density plot viewer to \"false\"", () => setProperty(page, "showYSelector", el("density plot viewer"), "false"));
      await session.step(76, "Then the \"y selector shown\" reading of density plot viewer should be \"false\"", () => readingReads(page, "y selector shown", el("density plot viewer"), "false"));
      await session.step(77, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["showXSelector","true"],["showYSelector","true"]]));
      await session.step(80, "Then the \"x selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "x selector shown", el("density plot viewer"), "true"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The description sits above the bins and Description Position moves it below them", async () => {
      await session.step(84, "When user sets \"description\" property of density plot viewer to \"Age against height\"", () => setProperty(page, "description", el("density plot viewer"), "Age against height"));
      await session.step(85, "Then the description of density plot viewer should be above its content", () => descriptionAbove(page, el("density plot viewer")));
      await session.step(86, "When user sets \"descriptionPosition\" property of density plot viewer to \"Bottom\"", () => setProperty(page, "descriptionPosition", el("density plot viewer"), "Bottom"));
      await session.step(87, "Then the description of density plot viewer should be below its content", () => descriptionBelow(page, el("density plot viewer")));
      await session.step(88, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["descriptionPosition","Top"],["description",""]]));
      await session.step(91, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own filter narrows what it bins and leaves the table alone", async () => {
      await session.step(94, "When user sets \"filter\" property of density plot viewer to \"${AGE} > 30\"", () => setProperty(page, "filter", el("density plot viewer"), "${AGE} > 30"));
      await session.step(95, "Then the \"rows shown\" reading of density plot viewer should be lower than before", () => readingLower(page, "rows shown", el("density plot viewer")));
      await session.step(96, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(97, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(98, "When user sets \"filter\" property of density plot viewer to \"\"", () => setProperty(page, "filter", el("density plot viewer"), ""));
      await session.step(99, "Then the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table moves what the plot bins", async () => {
      await session.step(103, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(104, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(105, "And the \"rows shown\" reading of density plot viewer should be lower than before", () => readingLower(page, "rows shown", el("density plot viewer")));
      await session.step(106, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(107, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(108, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(109, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(110, "And the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bound to another table, the plot bins that table's rows", async () => {
      await session.step(114, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(115, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(116, "When user sets \"table\" property of density plot viewer to \"spgi-100\"", () => setProperty(page, "table", el("density plot viewer"), "spgi-100"));
      await session.step(117, "Then density plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("density plot viewer"), "spgi-100"));
      await session.step(118, "And the \"rows shown\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("density plot viewer"), 1));
      await session.step(119, "And the \"rows shown\" reading of density plot viewer should differ from before", () => readingDiffers(page, "rows shown", el("density plot viewer")));
      await session.step(120, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
