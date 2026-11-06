/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/plan/plan-form.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {effectsArrived, figureHigher, figureLower, figureReads, figureSame, finishedRunning, rememberFigure} from '../../bindings/pmx.js';
import {openApp} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Plan's declaration form", () => {
  const session = feature(test, "features/plan/plan-form.feature", import.meta.url);
  test("Plan's declaration form", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens PMX Studio", () => openApp(page));
    await session.step(9, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(10, "And \"Experimental groups\" section should contain text \"Control\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Control"));
    await run.scenario("A loaded plan is not out of date and Run is on", async () => {
      await session.step(13, "Then stale note should be hidden", () => shouldBe(page, el("stale note"), "hidden"));
      await session.step(14, "And form run button should be enabled", () => shouldBe(page, el("form run button"), "enabled"));
      await session.step(15, "And \"Cell line\" input should have the value \"MC38\"", () => shouldHaveValue(page, el("\"Cell line\" input"), "MC38"));
      await session.step(16, "And \"Endpoint volume\" input should have the value \"1500\"", () => shouldHaveValue(page, el("\"Endpoint volume\" input"), "1500"));
      await session.step(17, "And \"Target power\" input should have the value \"80\"", () => shouldHaveValue(page, el("\"Target power\" input"), "80"));
    });
    await run.scenario("Cell line", async () => {
      await session.step(20, "When user selects \"CT26\" in \"Cell line\" input", () => selectIn(page, "CT26", el("\"Cell line\" input")));
      await session.step(21, "Then stale note should be visible", () => shouldBe(page, el("stale note"), "visible"));
      await session.step(22, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(23, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(24, "And stale note should be hidden", () => shouldBe(page, el("stale note"), "hidden"));
      await session.step(25, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(26, "Then \"Study protocol\" section should contain text \"CT26\"", () => shouldContainText(page, el("\"Study protocol\" section"), "CT26"));
      await session.step(27, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(28, "And user selects \"B16F10\" in \"Cell line\" input", () => selectIn(page, "B16F10", el("\"Cell line\" input")));
      await session.step(29, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(30, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(31, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(32, "Then \"Study protocol\" section should contain text \"B16F10\"", () => shouldContainText(page, el("\"Study protocol\" section"), "B16F10"));
      await session.step(33, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(34, "And user selects \"MC38\" in \"Cell line\" input", () => selectIn(page, "MC38", el("\"Cell line\" input")));
      await session.step(35, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(36, "Then the study should have finished running", () => finishedRunning(page));
    });
    await run.scenario("Endpoint volume", async () => {
      await session.step(39, "When user remembers the \"End of study\" figure", () => rememberFigure(page, "End of study"));
      await session.step(40, "And user enters \"800\" into \"Endpoint volume\" input", () => enterInto(page, "800", el("\"Endpoint volume\" input")));
      await session.step(41, "Then stale note should be visible", () => shouldBe(page, el("stale note"), "visible"));
      await session.step(42, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(43, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(44, "And the \"End of study\" figure should be lower than remembered", () => figureLower(page, "End of study"));
      await session.step(45, "When user enters \"100\" into \"Endpoint volume\" input", () => enterInto(page, "100", el("\"Endpoint volume\" input")));
      await session.step(46, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(47, "When user enters \"5000\" into \"Endpoint volume\" input", () => enterInto(page, "5000", el("\"Endpoint volume\" input")));
      await session.step(48, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(49, "When user clears \"Endpoint volume\" input", () => clearField(page, el("\"Endpoint volume\" input")));
      await session.step(50, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(51, "When user enters \"1500\" into \"Endpoint volume\" input", () => enterInto(page, "1500", el("\"Endpoint volume\" input")));
      await session.step(52, "Then form run button should be enabled", () => shouldBe(page, el("form run button"), "enabled"));
      await session.step(53, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(54, "Then the study should have finished running", () => finishedRunning(page));
    });
    await run.scenario("Target power", async () => {
      await session.step(57, "When user remembers the \"Analyzable\" figure", () => rememberFigure(page, "Analyzable"));
      await session.step(58, "And user enters \"90\" into \"Target power\" input", () => enterInto(page, "90", el("\"Target power\" input")));
      await session.step(59, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(60, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(61, "And the \"Analyzable\" figure should be higher than remembered", () => figureHigher(page, "Analyzable"));
      await session.step(62, "When user enters \"99\" into \"Target power\" input", () => enterInto(page, "99", el("\"Target power\" input")));
      await session.step(63, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(64, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(65, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
      await session.step(66, "And the \"Analyzable\" figure should read \"—\"", () => figureReads(page, "Analyzable", "—"));
      await session.step(67, "And \"Experimental groups\" section should contain text \"do not reach 99% power\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "do not reach 99% power"));
      await session.step(68, "When user clicks on \"Trade-offs\" tab", () => clickOn(page, el("\"Trade-offs\" tab")));
      await session.step(69, "Then \"Power by end of study\" section should contain text \"No group size\"", () => shouldContainText(page, el("\"Power by end of study\" section"), "No group size"));
      await session.step(70, "And \"Measurement schedule\" section should contain text \"No group size\"", () => shouldContainText(page, el("\"Measurement schedule\" section"), "No group size"));
      await session.step(71, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(72, "Then \"Study protocol\" section should contain text \"Not sized\"", () => shouldContainText(page, el("\"Study protocol\" section"), "Not sized"));
      await session.step(73, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(74, "When user enters \"10\" into \"Target power\" input", () => enterInto(page, "10", el("\"Target power\" input")));
      await session.step(75, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(76, "When user enters \"80\" into \"Target power\" input", () => enterInto(page, "80", el("\"Target power\" input")));
      await session.step(77, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(78, "Then the study should have finished running", () => finishedRunning(page));
    });
    await run.scenario("Animals per arm", async () => {
      await session.step(81, "When user enters \"12\" into \"Animals per arm\" input", () => enterInto(page, "12", el("\"Animals per arm\" input")));
      await session.step(82, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(83, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(84, "And the \"Analyzable\" figure should read \"24\"", () => figureReads(page, "Analyzable", "24"));
      await session.step(85, "And \"Experimental groups\" section should contain text \"Smallest detectable effect\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Smallest detectable effect"));
      await session.step(86, "And the smallest detectable effects should have arrived", () => effectsArrived(page));
      await session.step(87, "When user enters \"1\" into \"Animals per arm\" input", () => enterInto(page, "1", el("\"Animals per arm\" input")));
      await session.step(88, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(89, "When user clears \"Animals per arm\" input", () => clearField(page, el("\"Animals per arm\" input")));
      await session.step(90, "Then form run button should be enabled", () => shouldBe(page, el("form run button"), "enabled"));
      await session.step(91, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(92, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(93, "And the \"Analyzable\" figure should read \"16\"", () => figureReads(page, "Analyzable", "16"));
    });
    await run.scenario("End of study", async () => {
      await session.step(96, "When user enters \"25\" into \"End of study\" input", () => enterInto(page, "25", el("\"End of study\" input")));
      await session.step(97, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(98, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(99, "And the \"End of study\" figure should read \"Day 25\"", () => figureReads(page, "End of study", "Day 25"));
      await session.step(100, "When user enters \"3\" into \"End of study\" input", () => enterInto(page, "3", el("\"End of study\" input")));
      await session.step(101, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(102, "When user clears \"End of study\" input", () => clearField(page, el("\"End of study\" input")));
      await session.step(103, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(104, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(105, "And the \"End of study\" figure should read \"Day 18\"", () => figureReads(page, "End of study", "Day 18"));
    });
    await run.scenario("Advanced fields", async () => {
      await session.step(108, "Then \"Significance (α)\" input should be hidden", () => shouldBe(page, el("\"Significance (α)\" input"), "hidden"));
      await session.step(109, "When user clicks on advanced header", () => clickOn(page, el("advanced header")));
      await session.step(110, "Then \"Significance (α)\" input should be visible", () => shouldBe(page, el("\"Significance (α)\" input"), "visible"));
      await session.step(111, "And \"Lost before randomization\" input should be visible", () => shouldBe(page, el("\"Lost before randomization\" input"), "visible"));
      await session.step(112, "And \"Lost after randomization\" input should be visible", () => shouldBe(page, el("\"Lost after randomization\" input"), "visible"));
      await session.step(113, "And \"Extra animals for losses\" input should be visible", () => shouldBe(page, el("\"Extra animals for losses\" input"), "visible"));
      await session.step(114, "And \"Confidence\" input should be visible", () => shouldBe(page, el("\"Confidence\" input"), "visible"));
      await session.step(115, "And \"Volume at randomization\" input should be visible", () => shouldBe(page, el("\"Volume at randomization\" input"), "visible"));
      await session.step(116, "And \"Humane limit\" input should be visible", () => shouldBe(page, el("\"Humane limit\" input"), "visible"));
      await session.step(117, "And \"Measurement days\" input should be hidden", () => shouldBe(page, el("\"Measurement days\" input"), "hidden"));
      await session.step(118, "And \"Randomization day\" input should be hidden", () => shouldBe(page, el("\"Randomization day\" input"), "hidden"));
      await session.step(119, "And \"Compare measurement schedules\" input should be hidden", () => shouldBe(page, el("\"Compare measurement schedules\" input"), "hidden"));
    });
    await run.scenario("Losses", async () => {
      await session.step(122, "When user remembers the \"Animals to implant\" figure", () => rememberFigure(page, "Animals to implant"));
      await session.step(123, "And user remembers the \"To randomize\" figure", () => rememberFigure(page, "To randomize"));
      await session.step(124, "And user enters \"50\" into \"Lost before randomization\" input", () => enterInto(page, "50", el("\"Lost before randomization\" input")));
      await session.step(125, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(126, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(127, "And the \"Animals to implant\" figure should be higher than remembered", () => figureHigher(page, "Animals to implant"));
      await session.step(128, "And the \"To randomize\" figure should be as remembered", () => figureSame(page, "To randomize"));
      await session.step(129, "When user enters \"20\" into \"Lost before randomization\" input", () => enterInto(page, "20", el("\"Lost before randomization\" input")));
      await session.step(130, "And user enters \"30\" into \"Lost after randomization\" input", () => enterInto(page, "30", el("\"Lost after randomization\" input")));
      await session.step(131, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(132, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(133, "And the \"To randomize\" figure should be higher than remembered", () => figureHigher(page, "To randomize"));
      await session.step(134, "When user enters \"10\" into \"Lost after randomization\" input", () => enterInto(page, "10", el("\"Lost after randomization\" input")));
      await session.step(135, "And user selects \"Average loss only\" in \"Extra animals for losses\" input", () => selectIn(page, "Average loss only", el("\"Extra animals for losses\" input")));
      await session.step(136, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(137, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(138, "And the \"Animals to implant\" figure should be lower than remembered", () => figureLower(page, "Animals to implant"));
      await session.step(139, "When user selects \"Enough for every group\" in \"Extra animals for losses\" input", () => selectIn(page, "Enough for every group", el("\"Extra animals for losses\" input")));
      await session.step(140, "And user enters \"99\" into \"Confidence\" input", () => enterInto(page, "99", el("\"Confidence\" input")));
      await session.step(141, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(142, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(143, "And the \"Animals to implant\" figure should be higher than remembered", () => figureHigher(page, "Animals to implant"));
      await session.step(144, "When user enters \"90\" into \"Confidence\" input", () => enterInto(page, "90", el("\"Confidence\" input")));
      await session.step(145, "And user enters \"120\" into \"Lost before randomization\" input", () => enterInto(page, "120", el("\"Lost before randomization\" input")));
      await session.step(146, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(147, "When user enters \"20\" into \"Lost before randomization\" input", () => enterInto(page, "20", el("\"Lost before randomization\" input")));
      await session.step(148, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(149, "Then the study should have finished running", () => finishedRunning(page));
    });
    await run.scenario("Volume at randomization, significance and humane limit", async () => {
      await session.step(152, "When user remembers the \"End of study\" figure", () => rememberFigure(page, "End of study"));
      await session.step(153, "And user enters \"300\" into \"Volume at randomization\" input", () => enterInto(page, "300", el("\"Volume at randomization\" input")));
      await session.step(154, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(155, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(156, "And the \"End of study\" figure should be lower than remembered", () => figureLower(page, "End of study"));
      await session.step(157, "When user enters \"100\" into \"Volume at randomization\" input", () => enterInto(page, "100", el("\"Volume at randomization\" input")));
      await session.step(158, "And user enters \"0.01\" into \"Significance (α)\" input", () => enterInto(page, "0.01", el("\"Significance (α)\" input")));
      await session.step(159, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(160, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(161, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(162, "Then \"Study protocol\" section should contain text \"α = 0.01\"", () => shouldContainText(page, el("\"Study protocol\" section"), "α = 0.01"));
      await session.step(163, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(164, "And user enters \"0.05\" into \"Significance (α)\" input", () => enterInto(page, "0.05", el("\"Significance (α)\" input")));
      await session.step(165, "And user enters \"1000\" into \"Humane limit\" input", () => enterInto(page, "1000", el("\"Humane limit\" input")));
      await session.step(166, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(167, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(168, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(169, "Then \"Study protocol\" section should contain text \"1000 mm³\"", () => shouldContainText(page, el("\"Study protocol\" section"), "1000 mm³"));
      await session.step(170, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(171, "And user enters \"2000\" into \"Humane limit\" input", () => enterInto(page, "2000", el("\"Humane limit\" input")));
      await session.step(172, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(173, "Then the study should have finished running", () => finishedRunning(page));
    });
    run.finish();
  });
});
