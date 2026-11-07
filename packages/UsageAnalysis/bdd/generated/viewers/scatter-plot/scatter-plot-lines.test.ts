/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-lines.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#3695, core#1200]
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
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, hasArea, hasNoArea, lessInk, moreInk, noBalloons, noErrors, painted, propertyShouldBe, readingIs, repainted, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot trend lines", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-lines.feature", import.meta.url);
  test("Scatter plot trend lines", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#3695", "@edge-cases", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(23, "And the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
    await run.scenario("A regression line is drawn, on a logarithmic axis too", async () => {
      await session.step(26, "Then \"Show Regression Line\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Regression Line", el("scatter plot viewer"), "false"));
      await session.step(27, "And the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(28, "And scatter plot viewer should not have a \"regression line\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(29, "When user sets \"Show Regression Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Regression Line", el("scatter plot viewer"), "true"));
      await session.step(30, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(31, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(32, "And the \"regression line\" area of scatter plot viewer should be painted", () => areaPainted(page, "regression line", el("scatter plot viewer")));
      await session.step(33, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(34, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(35, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(37, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["Show Regression Line","false"]]), [["Y Axis Type","linear"],["Show Regression Line","false"]]);
      await session.step(40, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(41, "And scatter plot viewer should not have a \"regression line\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Regression Per Category fits one line per color, and the stats table is drawn, and one common line past 20 categories (public#3695)", async () => {
      await session.step(46, "Then \"Regression Per Category\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Regression Per Category", el("scatter plot viewer"), "true"));
      await session.step(47, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Show Regression Line","true"]]), [["Color","RACE"],["Show Regression Line","true"]]);
      await session.step(50, "Then the \"regression lines\" reading of scatter plot viewer should be 4", () => readingIs(page, "regression lines", el("scatter plot viewer"), 4));
      await session.step(51, "And \"Show Regression Line Equation\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Regression Line Equation", el("scatter plot viewer"), "true"));
      await session.step(52, "And \"Show Determination Coefficient\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Determination Coefficient", el("scatter plot viewer"), "true"));
      await session.step(53, "When user sets \"Regression Per Category\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Regression Per Category", el("scatter plot viewer"), "false"));
      await session.step(54, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(55, "And scatter plot viewer should have a \"regression stats\" area", () => hasArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(56, "And the \"regression stats\" area of scatter plot viewer should be painted", () => areaPainted(page, "regression stats", el("scatter plot viewer")));
      await session.step(57, "When user sets \"Show Regression Line Equation\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Show Regression Line Equation", el("scatter plot viewer"), "false"));
      await session.step(58, "Then scatter plot viewer should have a \"regression stats\" area", () => hasArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(59, "And the \"regression stats\" area of scatter plot viewer should be painted", () => areaPainted(page, "regression stats", el("scatter plot viewer")));
      await session.step(60, "When user sets \"Show Determination Coefficient\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Show Determination Coefficient", el("scatter plot viewer"), "false"));
      await session.step(61, "Then scatter plot viewer should not have a \"regression stats\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(62, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Regression Line Equation","true"],["Show Determination Coefficient","true"],["Regression Per Category","true"]]), [["Show Regression Line Equation","true"],["Show Determination Coefficient","true"],["Regression Per Category","true"]]);
      await session.step(66, "Then the \"regression lines\" reading of scatter plot viewer should be 4", () => readingIs(page, "regression lines", el("scatter plot viewer"), 4));
      await session.step(67, "When user sets \"Color\" property of scatter plot viewer to \"USUBJID\"", () => setProperty(page, "Color", el("scatter plot viewer"), "USUBJID"));
      await session.step(68, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(69, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(70, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Regression Line","false"],["Color",""]]), [["Show Regression Line","false"],["Color",""]]);
      await session.step(73, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A regression over a datetime axis takes a time unit", async () => {
      await session.step(77, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","STARTED"],["Show Regression Line","true"]]), [["X","STARTED"],["Show Regression Line","true"]]);
      await session.step(80, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(81, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(82, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(83, "When user sets \"X Map\" property of scatter plot viewer to \"year\"", () => setProperty(page, "X Map", el("scatter plot viewer"), "year"));
      await session.step(84, "Then \"X Map\" property of scatter plot viewer should be \"year\"", () => propertyShouldBe(page, "X Map", el("scatter plot viewer"), "year"));
      await session.step(85, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(87, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X Map",""],["X","WEIGHT"],["Show Regression Line","false"]]), [["X Map",""],["X","WEIGHT"],["Show Regression Line","false"]]);
      await session.step(91, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(92, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The moving average line, its window, its split and its deviation band", async () => {
      await session.step(96, "Then scatter plot viewer should not have a \"moving average\" area", () => hasNoArea(page, el("scatter plot viewer"), "moving average"));
      await session.step(97, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Moving Average Per Category","false"]]), [["Color","RACE"],["Moving Average Per Category","false"]]);
      await session.step(100, "And user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(101, "Then scatter plot viewer should have a \"moving average\" area", () => hasArea(page, el("scatter plot viewer"), "moving average"));
      await session.step(102, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(103, "When user sets \"Moving Average Window\" property of scatter plot viewer to \"200\"", () => setProperty(page, "Moving Average Window", el("scatter plot viewer"), "200"));
      await session.step(104, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(105, "When user sets \"Moving Average Per Category\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Moving Average Per Category", el("scatter plot viewer"), "true"));
      await session.step(106, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(107, "When user sets \"Show Moving Average Deviation\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Deviation", el("scatter plot viewer"), "true"));
      await session.step(108, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(109, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Moving Average Deviation","false"],["Moving Average Window","10"],["Show Moving Average Line","false"],["Color",""]]), [["Show Moving Average Deviation","false"],["Moving Average Window","10"],["Show Moving Average Line","false"],["Color",""]]);
      await session.step(114, "Then scatter plot viewer should not have a \"moving average\" area", () => hasNoArea(page, el("scatter plot viewer"), "moving average"));
      await session.step(115, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(116, "And the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Moving Average Window of 1 typed in the context panel is clamped to 2 (core#1200) [typed=1, clamped=2]", async () => {
      await session.step(121, "Given user sets \"Moving Average Window\" property of scatter plot viewer to \"10\"", () => setProperty(page, "Moving Average Window", el("scatter plot viewer"), "10"));
      await session.step(122, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(123, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(124, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(125, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(126, "When user enters \"1\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "1", el("\"Moving Average Window\" property in context panel")));
      await session.step(127, "Then \"Moving Average Window\" property of scatter plot viewer should be \"2\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "2"));
    });
    await run.scenario("A Moving Average Window of 101 typed in the context panel is clamped to 100 (core#1200) [typed=101, clamped=100]", async () => {
      await session.step(121, "Given user sets \"Moving Average Window\" property of scatter plot viewer to \"10\"", () => setProperty(page, "Moving Average Window", el("scatter plot viewer"), "10"));
      await session.step(122, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(123, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(124, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(125, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(126, "When user enters \"101\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "101", el("\"Moving Average Window\" property in context panel")));
      await session.step(127, "Then \"Moving Average Window\" property of scatter plot viewer should be \"100\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "100"));
    });
    run.finish();
  });
});
