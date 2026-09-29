/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-overlays.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaAtLeastTall, clickArea, hasArea, hasNoArea, lessInk, loadLayout, moreInk, noErrors, propertyShouldBe, readingIs, readingReads, repainted, reportsNoError, saveLayoutToServer, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart regression, moving average and formula lines", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-overlays.feature", import.meta.url);
  test("Line chart regression, moving average and formula lines", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(24, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(25, "And the \"regression lines\" reading of line chart viewer should be 0", () => readingIs(page, "regression lines", el("line chart viewer"), 0));
    await session.step(26, "And the \"moving average lines\" reading of line chart viewer should be 0", () => readingIs(page, "moving average lines", el("line chart viewer"), 0));
    await session.step(27, "And the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
    await session.step(28, "And line chart viewer should not have a \"regression line\" area", () => hasNoArea(page, el("line chart viewer"), "regression line"));
    await session.step(29, "And line chart viewer should not have a \"moving average\" area", () => hasNoArea(page, el("line chart viewer"), "moving average"));
    await session.step(30, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("The regression line is one fit over the whole series", async () => {
      await session.step(33, "When user sets \"showRegressionLine\" property of line chart viewer to \"true\"", () => setProperty(page, "showRegressionLine", el("line chart viewer"), "true"));
      await session.step(34, "Then the \"regression lines\" reading of line chart viewer should be 1", () => readingIs(page, "regression lines", el("line chart viewer"), 1));
      await session.step(35, "And line chart viewer should have a \"regression line\" area", () => hasArea(page, el("line chart viewer"), "regression line"));
      await session.step(36, "And line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(37, "When user sets \"showRegressionLine\" property of line chart viewer to \"false\"", () => setProperty(page, "showRegressionLine", el("line chart viewer"), "false"));
      await session.step(38, "Then the \"regression lines\" reading of line chart viewer should be 0", () => readingIs(page, "regression lines", el("line chart viewer"), 0));
      await session.step(39, "And line chart viewer should not have a \"regression line\" area", () => hasNoArea(page, el("line chart viewer"), "regression line"));
      await session.step(40, "And line chart viewer should have less ink than before", () => lessInk(page, el("line chart viewer")));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The R key toggles the regression line the chart says it binds", async () => {
      await session.step(44, "When user clicks on the \"plot\" area of line chart viewer", () => clickArea(page, "plot", el("line chart viewer")));
      await session.step(45, "And user presses R in line chart viewer", () => pressKeyIn(page, "R", el("line chart viewer")));
      await session.step(46, "Then the \"regression lines\" reading of line chart viewer should be 1", () => readingIs(page, "regression lines", el("line chart viewer"), 1));
      await session.step(47, "And line chart viewer should have a \"regression line\" area", () => hasArea(page, el("line chart viewer"), "regression line"));
      await session.step(48, "When user presses R in line chart viewer", () => pressKeyIn(page, "R", el("line chart viewer")));
      await session.step(49, "Then the \"regression lines\" reading of line chart viewer should be 0", () => readingIs(page, "regression lines", el("line chart viewer"), 0));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An R typed into the chart's own column selector is a letter, not the shortcut", async () => {
      await session.step(53, "When user picks \"NIBR logP\" in the \"x\" column selector of line chart viewer", () => pickInColumnSelector(page, "NIBR logP", "x", el("line chart viewer")));
      await session.step(54, "Then \"showRegressionLine\" property of line chart viewer should be \"false\"", () => propertyShouldBe(page, "showRegressionLine", el("line chart viewer"), "false"));
      await session.step(55, "And the \"regression lines\" reading of line chart viewer should be 0", () => readingIs(page, "regression lines", el("line chart viewer"), 0));
      await session.step(56, "When user picks \"CAST Idea ID\" in the \"x\" column selector of line chart viewer", () => pickInColumnSelector(page, "CAST Idea ID", "x", el("line chart viewer")));
      await session.step(57, "Then the \"x column\" reading of line chart viewer should be \"CAST Idea ID\"", () => readingReads(page, "x column", el("line chart viewer"), "CAST Idea ID"));
    });
    await run.scenario("The moving average is a series of its own, and the deviation band widens it", async () => {
      await session.step(60, "When user sets \"showMovingAverageLine\" property of line chart viewer to \"true\"", () => setProperty(page, "showMovingAverageLine", el("line chart viewer"), "true"));
      await session.step(61, "Then the \"moving average lines\" reading of line chart viewer should be 1", () => readingIs(page, "moving average lines", el("line chart viewer"), 1));
      await session.step(62, "And the \"moving average deviation\" reading of line chart viewer should be \"false\"", () => readingReads(page, "moving average deviation", el("line chart viewer"), "false"));
      await session.step(63, "And line chart viewer should have a \"moving average\" area", () => hasArea(page, el("line chart viewer"), "moving average"));
      await session.step(64, "And line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(65, "When user sets \"showMovingAverageDeviation\" property of line chart viewer to \"true\"", () => setProperty(page, "showMovingAverageDeviation", el("line chart viewer"), "true"));
      await session.step(66, "Then the \"moving average deviation\" reading of line chart viewer should be \"true\"", () => readingReads(page, "moving average deviation", el("line chart viewer"), "true"));
      await session.step(67, "And line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(68, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showMovingAverageDeviation","false"],["showMovingAverageLine","false"]]), [["showMovingAverageDeviation","false"],["showMovingAverageLine","false"]]);
      await session.step(71, "Then the \"moving average lines\" reading of line chart viewer should be 0", () => readingIs(page, "moving average lines", el("line chart viewer"), 0));
      await session.step(72, "And line chart viewer should not have a \"moving average\" area", () => hasNoArea(page, el("line chart viewer"), "moving average"));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split gives every category its own fit and its own moving average", async () => {
      await session.step(76, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showRegressionLine","true"],["showMovingAverageLine","true"]]), [["showRegressionLine","true"],["showMovingAverageLine","true"]]);
      await session.step(79, "Then the \"regression lines\" reading of line chart viewer should be 1", () => readingIs(page, "regression lines", el("line chart viewer"), 1));
      await session.step(80, "And the \"moving average lines\" reading of line chart viewer should be 1", () => readingIs(page, "moving average lines", el("line chart viewer"), 1));
      await session.step(81, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(82, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(83, "And the \"regression lines\" reading of line chart viewer should be 5", () => readingIs(page, "regression lines", el("line chart viewer"), 5));
      await session.step(84, "And the \"moving average lines\" reading of line chart viewer should be 5", () => readingIs(page, "moving average lines", el("line chart viewer"), 5));
      await session.step(85, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(86, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(87, "Then the \"regression lines\" reading of line chart viewer should be 1", () => readingIs(page, "regression lines", el("line chart viewer"), 1));
      await session.step(88, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showRegressionLine","false"],["showMovingAverageLine","false"]]), [["showRegressionLine","false"],["showMovingAverageLine","false"]]);
      await session.step(91, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A formula line and a formula band each get a hit area named after their title", async () => {
      await session.step(94, "When user sets \"formulaLines\" property of line chart viewer to '[{\"type\": \"line\", \"formula\": \"${Chemical Space X} = 5\", \"title\": \"const-line\", \"color\": \"#FF0000\"}, {\"type\": \"band\", \"formula\": \"${Chemical Space X} in(2, 8)\", \"title\": \"const-band\", \"color\": \"#00FF00\"}]'", () => setProperty(page, "formulaLines", el("line chart viewer"), "[{\"type\": \"line\", \"formula\": \"${Chemical Space X} = 5\", \"title\": \"const-line\", \"color\": \"#FF0000\"}, {\"type\": \"band\", \"formula\": \"${Chemical Space X} in(2, 8)\", \"title\": \"const-band\", \"color\": \"#00FF00\"}]"));
      await session.step(95, "Then the \"formula lines\" reading of line chart viewer should be 2", () => readingIs(page, "formula lines", el("line chart viewer"), 2));
      await session.step(96, "And line chart viewer should have a 'formula line \"const-line\"' area", () => hasArea(page, el("line chart viewer"), "formula line \"const-line\""));
      await session.step(97, "And line chart viewer should have a 'formula band \"const-band\"' area", () => hasArea(page, el("line chart viewer"), "formula band \"const-band\""));
      await session.step(98, "And line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(99, "And the 'formula band \"const-band\"' area of line chart viewer should be at least 10 pixels tall", () => areaAtLeastTall(page, "formula band \"const-band\"", el("line chart viewer"), 10));
      await session.step(100, "When user sets \"formulaLines\" property of line chart viewer to \"\"", () => setProperty(page, "formulaLines", el("line chart viewer"), ""));
      await session.step(101, "Then the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
      await session.step(102, "And line chart viewer should not have a 'formula line \"const-line\"' area", () => hasNoArea(page, el("line chart viewer"), "formula line \"const-line\""));
      await session.step(103, "And line chart viewer should not have a 'formula band \"const-band\"' area", () => hasNoArea(page, el("line chart viewer"), "formula band \"const-band\""));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The formula lines come back from a saved layout", async () => {
      await session.step(107, "When user sets \"formulaLines\" property of line chart viewer to '[{\"type\": \"line\", \"formula\": \"${Chemical Space X} = 5\", \"title\": \"const-line\", \"color\": \"#FF0000\"}, {\"type\": \"band\", \"formula\": \"${Chemical Space X} in(2, 8)\", \"title\": \"const-band\", \"color\": \"#00FF00\"}]'", () => setProperty(page, "formulaLines", el("line chart viewer"), "[{\"type\": \"line\", \"formula\": \"${Chemical Space X} = 5\", \"title\": \"const-line\", \"color\": \"#FF0000\"}, {\"type\": \"band\", \"formula\": \"${Chemical Space X} in(2, 8)\", \"title\": \"const-band\", \"color\": \"#00FF00\"}]"));
      await session.step(108, "Then the \"formula lines\" reading of line chart viewer should be 2", () => readingIs(page, "formula lines", el("line chart viewer"), 2));
      await session.step(109, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(110, "And user sets \"formulaLines\" property of line chart viewer to \"\"", () => setProperty(page, "formulaLines", el("line chart viewer"), ""));
      await session.step(111, "Then the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
      await session.step(112, "When user loads the saved layout", () => loadLayout(page));
      await session.step(113, "Then the \"formula lines\" reading of line chart viewer should be 2", () => readingIs(page, "formula lines", el("line chart viewer"), 2));
      await session.step(114, "And line chart viewer should have a 'formula line \"const-line\"' area", () => hasArea(page, el("line chart viewer"), "formula line \"const-line\""));
      await session.step(115, "And line chart viewer should have a 'formula band \"const-band\"' area", () => hasArea(page, el("line chart viewer"), "formula band \"const-band\""));
      await session.step(116, "When user sets \"formulaLines\" property of line chart viewer to \"\"", () => setProperty(page, "formulaLines", el("line chart viewer"), ""));
      await session.step(117, "Then the \"formula lines\" reading of line chart viewer should be 0", () => readingIs(page, "formula lines", el("line chart viewer"), 0));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
