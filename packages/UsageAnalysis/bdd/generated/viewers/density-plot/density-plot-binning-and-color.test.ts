/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/density-plot/density-plot-binning-and-color.feature
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
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaRepainted, noErrors, painted, readingAsRemembered, readingAtLeast, readingHigher, readingIs, readingLower, readingReads, readingSame, readingsEqual, rememberReading, repainted, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, typeInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Density plot binning, bin shape and the colour scale", () => {
  const session = feature(test, "features/viewers/density-plot/density-plot-binning-and-color.feature", import.meta.url);
  test("Density plot binning, bin shape and the colour scale", {tag: ["@journey", "@viewers", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"]]));
    await session.step(19, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
    await session.step(20, "And the \"bins\" reading of density plot viewer should be 50", () => readingIs(page, "bins", el("density plot viewer"), 50));
    await session.step(21, "And the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
    await session.step(22, "And the \"color transform\" reading of density plot viewer should be \"linear\"", () => readingReads(page, "color transform", el("density plot viewer"), "linear"));
    await session.step(23, "And the \"color scheme inverted\" reading of density plot viewer should be \"false\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "false"));
    await run.scenario("Fewer bins means fewer, fuller bins; more bins means more, emptier ones", async () => {
      await session.step(26, "Then the \"bins drawn\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(27, "When user remembers the \"bins drawn\" reading of density plot viewer", () => rememberReading(page, "bins drawn", el("density plot viewer")));
      await session.step(28, "And user sets \"bins\" property of density plot viewer to \"5\"", () => setProperty(page, "bins", el("density plot viewer"), "5"));
      await session.step(29, "Then the \"bins\" reading of density plot viewer should be 5", () => readingIs(page, "bins", el("density plot viewer"), 5));
      await session.step(30, "And the \"bins drawn\" reading of density plot viewer should be lower than before", () => readingLower(page, "bins drawn", el("density plot viewer")));
      await session.step(31, "And the \"max bin count\" reading of density plot viewer should be higher than before", () => readingHigher(page, "max bin count", el("density plot viewer")));
      await session.step(32, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(33, "When user sets \"bins\" property of density plot viewer to \"200\"", () => setProperty(page, "bins", el("density plot viewer"), "200"));
      await session.step(34, "Then the \"bins drawn\" reading of density plot viewer should be higher than before", () => readingHigher(page, "bins drawn", el("density plot viewer")));
      await session.step(35, "And the \"max bin count\" reading of density plot viewer should be lower than before", () => readingLower(page, "max bin count", el("density plot viewer")));
      await session.step(36, "When user sets \"bins\" property of density plot viewer to \"50\"", () => setProperty(page, "bins", el("density plot viewer"), "50"));
      await session.step(37, "Then the \"bins drawn\" reading of density plot viewer should be as remembered", () => readingAsRemembered(page, "bins drawn", el("density plot viewer")));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A single rectangular bin holds every row that was binned", async () => {
      await session.step(41, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["binShape","rectangle"],["bins","1"]]));
      await session.step(44, "Then the \"bins drawn\" reading of density plot viewer should be 1", () => readingIs(page, "bins drawn", el("density plot viewer"), 1));
      await session.step(45, "And the \"rows in densest bin\" and \"rows shown\" readings of density plot viewer should be the same", () => readingsEqual(page, "rows in densest bin", "rows shown", el("density plot viewer")));
      await session.step(46, "When user sets properties of density plot viewer:", () => setProperties(page, el("density plot viewer"), [["bins","50"],["binShape","hexagon"]]));
      await session.step(49, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The colour scale spans zero to the fullest bin, and the transform is reported", async () => {
      await session.step(52, "Then the \"color scale min\" reading of density plot viewer should be 0", () => readingIs(page, "color scale min", el("density plot viewer"), 0));
      await session.step(53, "And the \"color scale max\" and \"max bin count\" readings of density plot viewer should be the same", () => readingsEqual(page, "color scale max", "max bin count", el("density plot viewer")));
      await session.step(54, "When user sets \"colorTransformType\" property of density plot viewer to \"logarithmic\"", () => setProperty(page, "colorTransformType", el("density plot viewer"), "logarithmic"));
      await session.step(55, "Then the \"color transform\" reading of density plot viewer should be \"logarithmic\"", () => readingReads(page, "color transform", el("density plot viewer"), "logarithmic"));
      await session.step(56, "And the \"color scale min\" reading of density plot viewer should be 1", () => readingIs(page, "color scale min", el("density plot viewer"), 1));
      await session.step(57, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(58, "When user sets \"colorTransformType\" property of density plot viewer to \"linear\"", () => setProperty(page, "colorTransformType", el("density plot viewer"), "linear"));
      await session.step(59, "Then the \"color scale min\" reading of density plot viewer should be 0", () => readingIs(page, "color scale min", el("density plot viewer"), 0));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting the colour scheme repaints the bins and the scale", async () => {
      await session.step(63, "Then the \"color scheme inverted\" reading of density plot viewer should be \"false\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "false"));
      await session.step(64, "When user sets \"invertColorScheme\" property of density plot viewer to \"true\"", () => setProperty(page, "invertColorScheme", el("density plot viewer"), "true"));
      await session.step(65, "Then the \"color scheme inverted\" reading of density plot viewer should be \"true\"", () => readingReads(page, "color scheme inverted", el("density plot viewer"), "true"));
      await session.step(66, "And the \"color scale\" area of density plot viewer should have repainted", () => areaRepainted(page, "color scale", el("density plot viewer")));
      await session.step(67, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(68, "And the \"rows in densest bin\" reading of density plot viewer should be the same as before", () => readingSame(page, "rows in densest bin", el("density plot viewer")));
      await session.step(69, "When user sets \"invertColorScheme\" property of density plot viewer to \"false\"", () => setProperty(page, "invertColorScheme", el("density plot viewer"), "false"));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rectangles and hexagons bin the same rows into different shapes", async () => {
      await session.step(73, "Then the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
      await session.step(74, "When user remembers the \"rows shown\" reading of density plot viewer", () => rememberReading(page, "rows shown", el("density plot viewer")));
      await session.step(75, "And user sets \"binShape\" property of density plot viewer to \"rectangle\"", () => setProperty(page, "binShape", el("density plot viewer"), "rectangle"));
      await session.step(76, "Then the \"bin shape\" reading of density plot viewer should be \"rectangle\"", () => readingReads(page, "bin shape", el("density plot viewer"), "rectangle"));
      await session.step(77, "And the \"rows shown\" reading of density plot viewer should be as remembered", () => readingAsRemembered(page, "rows shown", el("density plot viewer")));
      await session.step(78, "And density plot viewer should have repainted", () => repainted(page, el("density plot viewer")));
      await session.step(79, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(80, "When user sets \"binShape\" property of density plot viewer to \"hexagon\"", () => setProperty(page, "binShape", el("density plot viewer"), "hexagon"));
      await session.step(81, "Then the \"bin shape\" reading of density plot viewer should be \"hexagon\"", () => readingReads(page, "bin shape", el("density plot viewer"), "hexagon"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The columns can be re-picked on the plot itself", async () => {
      await session.step(85, "When user picks \"HEIGHT\" in the \"y\" column selector of density plot viewer", () => pickInColumnSelector(page, "HEIGHT", "y", el("density plot viewer")));
      await session.step(86, "Then the \"y column\" reading of density plot viewer should be \"HEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "HEIGHT"));
      await session.step(87, "And the \"rows shown\" reading of density plot viewer should be 872", () => readingIs(page, "rows shown", el("density plot viewer"), 872));
      await session.step(88, "When user picks \"WEIGHT\" in the \"y\" column selector of density plot viewer", () => pickInColumnSelector(page, "WEIGHT", "y", el("density plot viewer")));
      await session.step(89, "Then the \"y column\" reading of density plot viewer should be \"WEIGHT\"", () => readingReads(page, "y column", el("density plot viewer"), "WEIGHT"));
      await session.step(90, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column selector offers numerical columns only (GROK-17118)", async () => {
      await session.step(94, "Then the \"x column\" reading of density plot viewer should be \"AGE\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE"));
      await session.step(95, "When user types \"SEX\" into the \"x\" column selector of density plot viewer", () => typeInColumnSelector(page, "SEX", "x", el("density plot viewer")));
      await session.step(96, "Then the \"x column\" reading of density plot viewer should be \"AGE\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE"));
      await session.step(97, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A calculated column can be binned like any other", async () => {
      await session.step(101, "Given user adds a calculated column \"AGE_PLUS\" with formula \"${AGE} + 1\"", () => addCalculated(page, "AGE_PLUS", "${AGE} + 1"));
      await session.step(102, "When user sets \"xColumnName\" property of density plot viewer to \"AGE_PLUS\"", () => setProperty(page, "xColumnName", el("density plot viewer"), "AGE_PLUS"));
      await session.step(103, "Then the \"x column\" reading of density plot viewer should be \"AGE_PLUS\"", () => readingReads(page, "x column", el("density plot viewer"), "AGE_PLUS"));
      await session.step(104, "And the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
      await session.step(105, "And density plot viewer should be painted", () => painted(page, el("density plot viewer")));
      await session.step(106, "When user sets \"xColumnName\" property of density plot viewer to \"AGE\"", () => setProperty(page, "xColumnName", el("density plot viewer"), "AGE"));
      await session.step(107, "And user removes \"AGE_PLUS\" column", () => removeColumn(page, "AGE_PLUS"));
      await session.step(108, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
