/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/functions/function-view-analyses.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.functions]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {analysisSwitch, runsDeleted} from '../../bindings/function-view.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldContainText, shouldHaveValue, shouldHaveValueBetween} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {packageInstalled, scriptOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The function view opens the analyses a script enables, and an accepted fit fills the form", () => {
  const session = feature(test, "features/functions/function-view-analyses.feature", import.meta.url);
  test("The function view opens the analyses a script enables, and an accepted fit fills the form", {tag: ["@journey", "@realizes:views.functions"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(16, "And the runs of the \"BDDFuncAnalyses{time}\" function are deleted at the end", () => runsDeleted(page, session.text("BDDFuncAnalyses{time}")));
    await session.step(17, "And a script \"BDDFuncAnalyses{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFuncAnalyses{time}"), "//language: javascript\n//input: int rows = 5 {caption: Rows}\n//input: double slope = 2 {caption: Slope}\n//output: double total {caption: Total}\n//meta.features: {\"fitting\": true, \"sens-analysis\": true}\ntotal = slope * rows * (rows - 1) / 2;"));
    await run.scenario("The ribbon offers the analyses the script enables", async () => {
      await session.step(28, "When user opens the address \"/func/Admin.BDDFuncAnalyses{time}?rows=5\"", () => openAddress(page, session.text("/func/Admin.BDDFuncAnalyses{time}?rows=5")));
      await session.step(29, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(30, "Then \"Total\" input should have value \"20\"", () => shouldHaveValue(page, el("\"Total\" input"), "20"));
      await session.step(31, "And \"Fitting\" icon should be visible", () => shouldBe(page, el("\"Fitting\" icon"), "visible"));
      await session.step(32, "And \"Sensitivity analysis\" icon should be visible", () => shouldBe(page, el("\"Sensitivity analysis\" icon"), "visible"));
    });
    await run.scenario("Sensitivity analysis runs the function over the range of an input", async () => {
      await session.step(35, "When user clicks on \"Sensitivity analysis\" icon", () => clickOn(page, el("\"Sensitivity analysis\" icon")));
      await session.step(36, "Then \"Samples\" input should be visible", () => shouldBe(page, el("\"Samples\" input"), "visible"));
      await session.step(37, "When user switches on \"Slope\" in the analysis form", () => analysisSwitch(page, "Slope"));
      await session.step(38, "And user enters \"1\" into \"Slope min\" input", () => enterInto(page, "1", el("\"Slope min\" input")));
      await session.step(39, "And user enters \"5\" into \"Slope max\" input", () => enterInto(page, "5", el("\"Slope max\" input")));
      await session.step(40, "And user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
      await session.step(41, "Then status bar should contain text \"Rows: 10\"", () => shouldContainText(page, el("status bar"), "Rows: 10"));
    });
    await run.scenario("Accepting a fit puts the fitted input into the function view's form", async () => {
      await session.step(44, "When user clicks on \"BDDFuncAnalyses{time}\" view tab", () => clickOn(page, el(session.text("\"BDDFuncAnalyses{time}\" view tab"))));
      await session.step(45, "And user clicks on \"Fitting\" icon", () => clickOn(page, el("\"Fitting\" icon")));
      await session.step(46, "Then \"Apply fitted parameters to the model\" icon should be visible", () => shouldBe(page, el("\"Apply fitted parameters to the model\" icon"), "visible"));
      await session.step(47, "When user switches on \"Slope\" in the analysis form", () => analysisSwitch(page, "Slope"));
      await session.step(48, "And user enters \"1\" into \"Slope (min)\" input", () => enterInto(page, "1", el("\"Slope (min)\" input")));
      await session.step(49, "And user enters \"20\" into \"Slope (max)\" input", () => enterInto(page, "20", el("\"Slope (max)\" input")));
      await session.step(50, "And user switches on \"Total\" in the analysis form", () => analysisSwitch(page, "Total"));
      await session.step(51, "And user enters \"100\" into \"Total\" input", () => enterInto(page, "100", el("\"Total\" input")));
      await session.step(52, "And user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
      await session.step(53, "Then status bar should contain text \"Rows: 1\"", () => shouldContainText(page, el("status bar"), "Rows: 1"));
      await session.step(54, "When user clicks on \"Apply fitted parameters to the model\" icon", () => clickOn(page, el("\"Apply fitted parameters to the model\" icon")));
      await session.step(55, "Then \"Select fitted parameters\" dialog should be visible", () => shouldBe(page, el("\"Select fitted parameters\" dialog"), "visible"));
      await session.step(56, "When user clicks on OK button in \"Select fitted parameters\" dialog", () => clickOn(page, el("OK button in \"Select fitted parameters\" dialog")));
      await session.step(57, "Then the \"BDDFuncAnalyses{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDFuncAnalyses{time}")));
      await session.step(58, "And \"Slope\" input should have a value between 9.9 and 10.1", () => shouldHaveValueBetween(page, el("\"Slope\" input"), 9.9, 10.1));
    });
    run.finish();
  });
});
