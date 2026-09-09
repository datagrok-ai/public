/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-lines.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, hasArea, hasNoArea, lessInk, moreInk, noBalloons, noErrors, painted, propertyShouldBe, readingIs, repainted, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot trend lines", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-lines.feature", import.meta.url);
  test("Scatter plot trend lines", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]));
    await session.step(18, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(19, "And the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await run.scenario("A regression line is drawn, on a logarithmic axis too", async () => {
      await session.step(22, "Then \"Show Regression Line\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Regression Line", el("scatter plot viewer"), "false"));
      await session.step(23, "And the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(24, "And scatter plot viewer should not have a \"regression line\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(25, "When user sets \"Show Regression Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Regression Line", el("scatter plot viewer"), "true"));
      await session.step(26, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(27, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(28, "And the \"regression line\" area of scatter plot viewer should be painted", () => areaPainted(page, "regression line", el("scatter plot viewer")));
      await session.step(29, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(30, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(31, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(33, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["Show Regression Line","false"]]));
      await session.step(36, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(37, "And scatter plot viewer should not have a \"regression line\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Regression Per Category fits one line per color, and the equation table is drawn", async () => {
      await session.step(41, "Then \"Regression Per Category\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Regression Per Category", el("scatter plot viewer"), "true"));
      await session.step(42, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Show Regression Line","true"]]));
      await session.step(45, "Then the \"regression lines\" reading of scatter plot viewer should be 4", () => readingIs(page, "regression lines", el("scatter plot viewer"), 4));
      await session.step(46, "And \"Show Regression Line Equation\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Regression Line Equation", el("scatter plot viewer"), "true"));
      await session.step(47, "When user sets \"Regression Per Category\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Regression Per Category", el("scatter plot viewer"), "false"));
      await session.step(48, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(49, "And scatter plot viewer should have a \"regression stats\" area", () => hasArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(50, "And the \"regression stats\" area of scatter plot viewer should be painted", () => areaPainted(page, "regression stats", el("scatter plot viewer")));
      await session.step(51, "When user sets \"Show Regression Line Equation\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Show Regression Line Equation", el("scatter plot viewer"), "false"));
      await session.step(52, "Then scatter plot viewer should not have a \"regression stats\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(53, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Regression Line Equation","true"],["Regression Per Category","true"]]));
      await session.step(56, "Then the \"regression lines\" reading of scatter plot viewer should be 4", () => readingIs(page, "regression lines", el("scatter plot viewer"), 4));
      await session.step(57, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Regression Line","false"],["Color",""]]));
      await session.step(60, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A regression over a datetime axis takes a time unit", async () => {
      await session.step(64, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","STARTED"],["Show Regression Line","true"]]));
      await session.step(67, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(68, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(70, "When user sets \"X Map\" property of scatter plot viewer to \"year\"", () => setProperty(page, "X Map", el("scatter plot viewer"), "year"));
      await session.step(71, "Then \"X Map\" property of scatter plot viewer should be \"year\"", () => propertyShouldBe(page, "X Map", el("scatter plot viewer"), "year"));
      await session.step(72, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X Map",""],["X","WEIGHT"],["Show Regression Line","false"]]));
      await session.step(78, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(79, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The moving average line, its window, its split and its deviation band", async () => {
      await session.step(83, "Then scatter plot viewer should not have a \"moving average\" area", () => hasNoArea(page, el("scatter plot viewer"), "moving average"));
      await session.step(84, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Moving Average Per Category","false"]]));
      await session.step(87, "And user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(88, "Then scatter plot viewer should have a \"moving average\" area", () => hasArea(page, el("scatter plot viewer"), "moving average"));
      await session.step(89, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(90, "When user sets \"Moving Average Window\" property of scatter plot viewer to \"200\"", () => setProperty(page, "Moving Average Window", el("scatter plot viewer"), "200"));
      await session.step(91, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(92, "When user sets \"Moving Average Per Category\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Moving Average Per Category", el("scatter plot viewer"), "true"));
      await session.step(93, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(94, "When user sets \"Show Moving Average Deviation\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Deviation", el("scatter plot viewer"), "true"));
      await session.step(95, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(96, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Moving Average Deviation","false"],["Moving Average Window","10"],["Show Moving Average Line","false"],["Color",""]]));
      await session.step(101, "Then scatter plot viewer should not have a \"moving average\" area", () => hasNoArea(page, el("scatter plot viewer"), "moving average"));
      await session.step(102, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(103, "And the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
