/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/density-plot/density-plot-binning-and-color.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.density-plot]
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
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, areaRepainted, loadLayout, noErrors, painted, readingAsRemembered, readingAtLeast, readingHigher, readingIs, readingLower, readingNotAsRemembered, readingReads, readingSame, readingsEqual, rememberReading, repainted, reportsNoError, saveLayout, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, typeInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Density plot binning, bin shape and the colour scale", () => {
  const session = feature(test, "features/viewers/density-plot/density-plot-binning-and-color.feature", import.meta.url);
  test("Density plot binning, bin shape and the colour scale", {tag: ["@journey", "@viewers", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(23, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
    await session.step(24, "And the \"bins\" reading of density plot viewer should be 50", () => readingIs(page, "bins", el("density plot viewer"), 50));
    await session.step(25, "And the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
    await session.step(26, "And the \"color transform\" reading of density plot viewer should be \"linear\"", () => readingReads(page, "color transform", el("density plot viewer"), "linear"));
    await session.step(27, "And the \"color scheme inverted\" reading of density plot viewer should be \"false\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "false"));
    await run.scenario("Fewer bins means fewer, fuller bins; more bins means more, emptier ones", async () => {
      await session.step(30, "Then the \"bins drawn\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(31, "When user remembers the \"bins drawn\" reading of density plot viewer", () => rememberReading(page, "bins drawn", el("density plot viewer")));
      await session.step(32, "And user sets \"bins\" property of density plot viewer to \"5\"", () => setProperty(page, "bins", el("density plot viewer"), "5"));
      await session.step(33, "Then the \"bins\" reading of density plot viewer should be 5", () => readingIs(page, "bins", el("density plot viewer"), 5));
      await session.step(34, "And the \"bins drawn\" reading of density plot viewer should be lower than before", () => readingLower(page, "bins drawn", el("density plot viewer")));
      await session.step(35, "And the \"max bin count\" reading of density plot viewer should be higher than before", () => readingHigher(page, "max bin count", el("density plot viewer")));
      await session.step(36, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(37, "When user sets \"bins\" property of density plot viewer to \"200\"", () => setProperty(page, "bins", el("density plot viewer"), "200"));
      await session.step(38, "Then the \"bins drawn\" reading of density plot viewer should be higher than before", () => readingHigher(page, "bins drawn", el("density plot viewer")));
      await session.step(39, "And the \"max bin count\" reading of density plot viewer should be lower than before", () => readingLower(page, "max bin count", el("density plot viewer")));
      await session.step(40, "When user sets \"bins\" property of density plot viewer to \"50\"", () => setProperty(page, "bins", el("density plot viewer"), "50"));
      await session.step(41, "Then the \"bins drawn\" reading of density plot viewer should be as remembered", () => readingAsRemembered(page, "bins drawn", el("density plot viewer")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A single rectangular bin holds every row that was binned", async () => {
      await session.step(45, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["binShape","rectangle"],["bins","1"]]), [["binShape","rectangle"],["bins","1"]]);
      await session.step(48, "Then the \"bins drawn\" reading of density plot viewer should be 1", () => readingIs(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(49, "And the \"rows in densest bin\" reading of density plot viewer should be 1000", () => readingIs(page, "rows in densest bin", el("density plot viewer"), 1000));
      await session.step(50, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["bins","50"],["binShape","hexagon"]]), [["bins","50"],["binShape","hexagon"]]);
      await session.step(53, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The colour scale spans zero to the fullest bin, and the transform is reported", async () => {
      await session.step(56, "Then the \"color scale min\" reading of density plot viewer should be 0", () => readingIs(page, "color scale min", el("density plot viewer"), 0));
      await session.step(57, "And the \"color scale max\" and \"max bin count\" readings of density plot viewer should be the same", () => readingsEqual(page, "color scale max", "max bin count", el("density plot viewer")));
      await session.step(58, "When user sets \"colorTransformType\" property of density plot viewer to \"logarithmic\"", () => setProperty(page, "colorTransformType", el("density plot viewer"), "logarithmic"));
      await session.step(59, "Then the \"color transform\" reading of density plot viewer should be \"logarithmic\"", () => readingReads(page, "color transform", el("density plot viewer"), "logarithmic"));
      await session.step(60, "And the \"color scale min\" reading of density plot viewer should be 1", () => readingIs(page, "color scale min", el("density plot viewer"), 1));
      await session.step(61, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(62, "When user sets \"colorTransformType\" property of density plot viewer to \"linear\"", () => setProperty(page, "colorTransformType", el("density plot viewer"), "linear"));
      await session.step(63, "Then the \"color scale min\" reading of density plot viewer should be 0", () => readingIs(page, "color scale min", el("density plot viewer"), 0));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting the colour scheme repaints the bins and the scale", async () => {
      await session.step(67, "Then the \"color scheme inverted\" reading of density plot viewer should be \"false\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "false"));
      await session.step(68, "When user sets \"invertColorScheme\" property of density plot viewer to \"true\"", () => setProperty(page, "invertColorScheme", el("density plot viewer"), "true"));
      await session.step(69, "Then the \"color scheme inverted\" reading of density plot viewer should be \"true\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "true"));
      await session.step(70, "And the \"color scale\" area of density plot viewer should have repainted", () => areaRepainted(page, "color scale", el("density plot viewer")));
      await session.step(71, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(72, "And the \"rows in densest bin\" reading of density plot viewer should be the same as before", () => readingSame(page, "rows in densest bin", el("density plot viewer")));
      await session.step(73, "When user sets \"invertColorScheme\" property of density plot viewer to \"false\"", () => setProperty(page, "invertColorScheme", el("density plot viewer"), "false"));
      await session.step(74, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rectangles and hexagons bin the same rows into different shapes", async () => {
      await session.step(77, "Then the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
      await session.step(78, "When user remembers the \"rows shown\" reading of density plot viewer", () => rememberReading(page, "rows shown", el("density plot viewer")));
      await session.step(79, "And user sets \"binShape\" property of density plot viewer to \"rectangle\"", () => setProperty(page, "binShape", el("density plot viewer"), "rectangle"));
      await session.step(80, "Then the \"bin shape\" reading of density plot viewer should be \"rectangle\"", () => readingReads(page, "bin shape", el("density plot viewer"), "rectangle"));
      await session.step(81, "And the \"rows shown\" reading of density plot viewer should be as remembered", () => readingAsRemembered(page, "rows shown", el("density plot viewer")));
      await session.step(82, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(83, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(84, "When user sets \"binShape\" property of density plot viewer to \"hexagon\"", () => setProperty(page, "binShape", el("density plot viewer"), "hexagon"));
      await session.step(85, "Then the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The columns can be re-picked on the plot itself", async () => {
      await session.step(89, "When user picks \"HEIGHT\" in the \"y\" column selector of density plot viewer", () => pickInColumnSelector(page, "HEIGHT", "y", el("density plot viewer")));
      await session.step(90, "Then the \"y column\" reading of density plot viewer should be \"HEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "HEIGHT"));
      await session.step(91, "And the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(92, "When user picks \"WEIGHT\" in the \"y\" column selector of density plot viewer", () => pickInColumnSelector(page, "WEIGHT", "y", el("density plot viewer")));
      await session.step(93, "Then the \"y column\" reading of density plot viewer should be \"WEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "WEIGHT"));
      await session.step(94, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column selector offers numerical columns only (GROK-17118)", async () => {
      await session.step(98, "Then the \"x column\" reading of density plot viewer should be \"AGE\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE"));
      await session.step(99, "When user types \"SEX\" into the \"x\" column selector of density plot viewer", () => typeInColumnSelector(page, "SEX", "x", el("density plot viewer")));
      await session.step(100, "Then the \"x column\" reading of density plot viewer should be \"AGE\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE"));
      await session.step(101, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A calculated column can be binned like any other", async () => {
      await session.step(105, "Given user adds a calculated column \"AGE_PLUS\" with formula \"${AGE} + 1\"", () => addCalculated(page, "AGE_PLUS", "${AGE} + 1"));
      await session.step(106, "When user sets \"xColumnName\" property of density plot viewer to \"AGE_PLUS\"", () => setProperty(page, "xColumnName", el("density plot viewer"), "AGE_PLUS"));
      await session.step(107, "Then the \"x column\" reading of density plot viewer should be \"AGE_PLUS\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE_PLUS"));
      await session.step(108, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(109, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(110, "When user sets \"xColumnName\" property of density plot viewer to \"AGE\"", () => setProperty(page, "xColumnName", el("density plot viewer"), "AGE"));
      await session.step(111, "And user removes \"AGE_PLUS\" column", () => removeColumn(page, "AGE_PLUS"));
      await session.step(112, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column with one value bins into a single strip in either shape", async () => {
      await session.step(115, "Given user adds a calculated column \"CONSTANT\" with formula \"${AGE} * 0 + 42\"", () => addCalculated(page, "CONSTANT", "${AGE} * 0 + 42"));
      await session.step(116, "When user picks \"CONSTANT\" in the \"x\" column selector of density plot viewer", () => pickInColumnSelector(page, "CONSTANT", "x", el("density plot viewer")));
      await session.step(117, "Then the \"x column\" reading of density plot viewer should be \"CONSTANT\"", () => readingReads(page, "x column", el("density plot viewer"), "CONSTANT"));
      await session.step(118, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(119, "And the \"bins drawn\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(120, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(121, "When user sets \"binShape\" property of density plot viewer to \"rectangle\"", () => setProperty(page, "binShape", el("density plot viewer"), "rectangle"));
      await session.step(122, "Then the \"bin shape\" reading of density plot viewer should be \"rectangle\"", () => readingReads(page, "bin shape", el("density plot viewer"), "rectangle"));
      await session.step(123, "And the \"bins drawn\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(124, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(125, "And density plot viewer should report no error", () => reportsNoError(page, el("density plot viewer")));
      await session.step(126, "When user sets \"binShape\" property of density plot viewer to \"hexagon\"", () => setProperty(page, "binShape", el("density plot viewer"), "hexagon"));
      await session.step(127, "Then the \"bins drawn\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(128, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(129, "When user picks \"AGE\" in the \"x\" column selector of density plot viewer", () => pickInColumnSelector(page, "AGE", "x", el("density plot viewer")));
      await session.step(130, "And user removes \"CONSTANT\" column", () => removeColumn(page, "CONSTANT"));
      await session.step(131, "Then the \"x column\" reading of density plot viewer should be \"AGE\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE"));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The plot comes back as it was drawn from a saved layout and from a project", async () => {
      await session.step(135, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["yColumnName","HEIGHT"],["bins","200"],["binShape","rectangle"],["invertColorScheme","true"]]), [["yColumnName","HEIGHT"],["bins","200"],["binShape","rectangle"],["invertColorScheme","true"]]);
      await session.step(140, "Then the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(141, "And the \"bins\" reading of density plot viewer should be 200", () => readingIs(page, "bins", el("density plot viewer"), 200));
      await session.step(142, "When user remembers the \"bins drawn\" reading of density plot viewer", () => rememberReading(page, "bins drawn", el("density plot viewer")));
      await session.step(143, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(144, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(145, "And user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["yColumnName","WEIGHT"],["bins","50"],["binShape","hexagon"],["invertColorScheme","false"]]), [["yColumnName","WEIGHT"],["bins","50"],["binShape","hexagon"],["invertColorScheme","false"]]);
      await session.step(150, "Then the \"bins drawn\" reading of density plot viewer should not be as remembered", () => readingNotAsRemembered(page, "bins drawn", el("density plot viewer")));
      await session.step(151, "When user loads the saved layout", () => loadLayout(page));
      await session.step(152, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(153, "And the \"y column\" reading of density plot viewer should be \"HEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "HEIGHT"));
      await session.step(154, "And the \"bins\" reading of density plot viewer should be 200", () => readingIs(page, "bins", el("density plot viewer"), 200));
      await session.step(155, "And the \"bin shape\" reading of density plot viewer should be \"rectangle\"", () => readingReads(page, "bin shape", el("density plot viewer"), "rectangle"));
      await session.step(156, "And the \"color scheme inverted\" reading of density plot viewer should be \"true\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "true"));
      await session.step(157, "And the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(158, "And the \"bins drawn\" reading of density plot viewer should be as remembered", () => readingAsRemembered(page, "bins drawn", el("density plot viewer")));
      await session.step(159, "When user saves the current view as project \"bdd-density-plot-peak\"", () => saveAsProject(page, "bdd-density-plot-peak"));
      await session.step(160, "And user closes all views", () => closeAllViews(page));
      await session.step(161, "And user opens the \"bdd-density-plot-peak\" project", () => openProject(page, "bdd-density-plot-peak"));
      await session.step(162, "Then the open tableview should have 1 density plot viewer", () => viewerCount(page, 1, "density plot"));
      await session.step(163, "And the \"y column\" reading of density plot viewer should be \"HEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "HEIGHT"));
      await session.step(164, "And the \"bin shape\" reading of density plot viewer should be \"rectangle\"", () => readingReads(page, "bin shape", el("density plot viewer"), "rectangle"));
      await session.step(165, "And the \"color scheme inverted\" reading of density plot viewer should be \"true\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "true"));
      await session.step(166, "And the \"bins drawn\" reading of density plot viewer should be as remembered", () => readingAsRemembered(page, "bins drawn", el("density plot viewer")));
      await session.step(167, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
