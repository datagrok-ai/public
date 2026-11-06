/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/analyze.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {figureLower, figureReads, finishedRunning, pickLastAnimal, rememberFigure, rememberTwin, twinDiffers} from '../../bindings/pmx.js';
import {openAnalyze} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, clipboardContains, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {infoBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Analyze", () => {
  const session = feature(test, "features/analyze/analyze.feature", import.meta.url);
  test("Analyze", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens PMX Studio on Analyze", () => openAnalyze(page));
    await run.scenario("Start page", async () => {
      await session.step(11, "Then analyze start should be visible", () => shouldBe(page, el("analyze start"), "visible"));
      await session.step(12, "And \"Analyze an example study\" button should be enabled", () => shouldBe(page, el("\"Analyze an example study\" button"), "enabled"));
      await session.step(13, "And form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
    });
    await run.scenario("Example study, Efficacy tab", async () => {
      await session.step(16, "When user clicks on \"Analyze an example study\" button", () => clickOn(page, el("\"Analyze an example study\" button")));
      await session.step(17, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(18, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
      await session.step(19, "And analyze start should be hidden", () => shouldBe(page, el("analyze start"), "hidden"));
      await session.step(20, "And stale note should be hidden", () => shouldBe(page, el("stale note"), "hidden"));
      await session.step(21, "And \"Experimental groups\" section should contain text \"Reference\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Reference"));
      await session.step(22, "And \"Experimental groups\" section should contain text \"TGI 70%\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "TGI 70%"));
      await session.step(23, "And \"Observed tumor growth\" section should be visible", () => shouldBe(page, el("\"Observed tumor growth\" section"), "visible"));
      await session.step(24, "When user clicks on linear axis choice", () => clickOn(page, el("linear axis choice")));
      await session.step(25, "Then visible volume axis should have the value \"Linear\"", () => shouldHaveValue(page, el("visible volume axis"), "Linear"));
    });
    await run.scenario("Each animal tab", async () => {
      await session.step(28, "When user clicks on \"Each animal\" tab", () => clickOn(page, el("\"Each animal\" tab")));
      await session.step(29, "Then \"Each animal against its untreated twin\" section should be visible", () => shouldBe(page, el("\"Each animal against its untreated twin\" section"), "visible"));
      await session.step(30, "And animal list should contain text \"respond\"", () => shouldContainText(page, el("animal list"), "respond"));
      await session.step(31, "And twin effect line should contain text \"than its untreated twin\"", () => shouldContainText(page, el("twin effect line"), "than its untreated twin"));
      await session.step(32, "When user remembers the twin line", () => rememberTwin(page));
      await session.step(33, "And user clicks on the last animal row", () => pickLastAnimal(page));
      await session.step(34, "Then the twin line should differ from remembered", () => twinDiffers(page));
    });
    await run.scenario("Tolerability tab", async () => {
      await session.step(37, "When user clicks on \"Tolerability\" tab", () => clickOn(page, el("\"Tolerability\" tab")));
      await session.step(38, "Then \"Weight loss\" section should be visible", () => shouldBe(page, el("\"Weight loss\" section"), "visible"));
      await session.step(39, "And \"Body weight change from day 0\" section should contain text \"as weighed\"", () => shouldContainText(page, el("\"Body weight change from day 0\" section"), "as weighed"));
      await session.step(40, "When user checks \"Tumor-corrected\" input", () => check(page, el("\"Tumor-corrected\" input")));
      await session.step(41, "Then \"Body weight change from day 0\" section should contain text \"without the tumor\"", () => shouldContainText(page, el("\"Body weight change from day 0\" section"), "without the tumor"));
      await session.step(42, "When user unchecks \"Tumor-corrected\" input", () => uncheck(page, el("\"Tumor-corrected\" input")));
      await session.step(43, "Then \"Body weight change from day 0\" section should contain text \"as weighed\"", () => shouldContainText(page, el("\"Body weight change from day 0\" section"), "as weighed"));
    });
    await run.scenario("Plan vs observed tab", async () => {
      await session.step(46, "When user clicks on \"Plan vs observed\" tab", () => clickOn(page, el("\"Plan vs observed\" tab")));
      await session.step(47, "Then \"Effect of each treatment\" section should be visible", () => shouldBe(page, el("\"Effect of each treatment\" section"), "visible"));
      await session.step(48, "And \"Control growth against the prediction\" section should be visible", () => shouldBe(page, el("\"Control growth against the prediction\" section"), "visible"));
      await session.step(49, "And \"End of the study\" section should contain text \"median control\"", () => shouldContainText(page, el("\"End of the study\" section"), "median control"));
      await session.step(50, "And \"Contribute to the atlas\" section should contain text \"Coming soon\"", () => shouldContainText(page, el("\"Contribute to the atlas\" section"), "Coming soon"));
    });
    await run.scenario("Statistics tab", async () => {
      await session.step(53, "When user clicks on \"Statistics\" tab", () => clickOn(page, el("\"Statistics\" tab")));
      await session.step(54, "Then \"Results summary\" section should contain text \"evaluated on day\"", () => shouldContainText(page, el("\"Results summary\" section"), "evaluated on day"));
      await session.step(55, "And \"Growth rates\" section should be visible", () => shouldBe(page, el("\"Growth rates\" section"), "visible"));
      await session.step(56, "And \"Control check\" section should contain text \"control measurements\"", () => shouldContainText(page, el("\"Control check\" section"), "control measurements"));
      await session.step(57, "And \"Analysis protocol\" section should contain text \"Primary comparison\"", () => shouldContainText(page, el("\"Analysis protocol\" section"), "Primary comparison"));
      await session.step(58, "When user clicks on first copy button", () => clickOn(page, el("first copy button")));
      await session.step(59, "Then an info balloon containing \"Results summary copied\" should have been shown", () => infoBalloonText(page, "Results summary copied"));
      await session.step(60, "And the clipboard should contain text \"evaluated on day\"", () => clipboardContains(page, "evaluated on day"));
    });
    await run.scenario("The form's fields mark the analysis out of date and Run brings them in", async () => {
      await session.step(63, "When user enters \"1000\" into \"Endpoint volume\" input", () => enterInto(page, "1000", el("\"Endpoint volume\" input")));
      await session.step(64, "Then stale note should be visible", () => shouldBe(page, el("stale note"), "visible"));
      await session.step(65, "When user remembers the \"Evaluation day\" figure", () => rememberFigure(page, "Evaluation day"));
      await session.step(66, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(67, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(68, "And stale note should be hidden", () => shouldBe(page, el("stale note"), "hidden"));
      await session.step(69, "And the \"Evaluation day\" figure should be lower than remembered", () => figureLower(page, "Evaluation day"));
      await session.step(70, "When user enters \"100\" into \"Endpoint volume\" input", () => enterInto(page, "100", el("\"Endpoint volume\" input")));
      await session.step(71, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(72, "When user clears \"Endpoint volume\" input", () => clearField(page, el("\"Endpoint volume\" input")));
      await session.step(73, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(74, "When user enters \"1500\" into \"Endpoint volume\" input", () => enterInto(page, "1500", el("\"Endpoint volume\" input")));
      await session.step(75, "And user selects \"CT26\" in \"Cell line\" input", () => selectIn(page, "CT26", el("\"Cell line\" input")));
      await session.step(76, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(77, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(78, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
      await session.step(79, "When user selects \"MC38\" in \"Cell line\" input", () => selectIn(page, "MC38", el("\"Cell line\" input")));
      await session.step(80, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(81, "Then the study should have finished running", () => finishedRunning(page));
    });
    await run.scenario("Advanced fields", async () => {
      await session.step(84, "When user clicks on advanced header", () => clickOn(page, el("advanced header")));
      await session.step(85, "Then \"Evaluation day\" input should be visible", () => shouldBe(page, el("\"Evaluation day\" input"), "visible"));
      await session.step(86, "And \"Significance (α)\" input should be visible", () => shouldBe(page, el("\"Significance (α)\" input"), "visible"));
      await session.step(87, "And \"Humane limit, weight loss\" input should be visible", () => shouldBe(page, el("\"Humane limit, weight loss\" input"), "visible"));
      await session.step(88, "And \"Humane limit, tumor\" input should be visible", () => shouldBe(page, el("\"Humane limit, tumor\" input"), "visible"));
      await session.step(89, "And \"Tumor density\" input should be visible", () => shouldBe(page, el("\"Tumor density\" input"), "visible"));
      await session.step(90, "When user enters \"9\" into \"Evaluation day\" input", () => enterInto(page, "9", el("\"Evaluation day\" input")));
      await session.step(91, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(92, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(93, "And the \"Evaluation day\" figure should read \"Day 9\"", () => figureReads(page, "Evaluation day", "Day 9"));
      await session.step(94, "When user clears \"Evaluation day\" input", () => clearField(page, el("\"Evaluation day\" input")));
      await session.step(95, "And user enters \"1.2\" into \"Tumor density\" input", () => enterInto(page, "1.2", el("\"Tumor density\" input")));
      await session.step(96, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(97, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(98, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
    });
    await run.scenario("An evaluation day that cannot be analyzed says why", async () => {
      await session.step(101, "When user enters \"5\" into \"Evaluation day\" input", () => enterInto(page, "5", el("\"Evaluation day\" input")));
      await session.step(102, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(103, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(104, "And study error should contain text \"choose a measurement day\"", () => shouldContainText(page, el("study error"), "choose a measurement day"));
      await session.step(105, "When user enters \"0\" into \"Evaluation day\" input", () => enterInto(page, "0", el("\"Evaluation day\" input")));
      await session.step(106, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(107, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(108, "And study error should contain text \"choose a measurement day\"", () => shouldContainText(page, el("study error"), "choose a measurement day"));
      await session.step(109, "When user clears \"Evaluation day\" input", () => clearField(page, el("\"Evaluation day\" input")));
      await session.step(110, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(111, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(112, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
    });
    run.finish();
  });
});
