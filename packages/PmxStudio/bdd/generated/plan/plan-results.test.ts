/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/plan/plan-results.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {finishedRunning, noSavedPlan, planIsSaved} from '../../bindings/pmx.js';
import {openApp} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardContains, downloadThrough, downloadedContains, enterInto, fileMatchingDownloaded, selectIn, shouldBe, shouldContainText, shouldHaveValue, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {infoBalloonText, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Plan's results", () => {
  const session = feature(test, "features/plan/plan-results.feature", import.meta.url);
  test("Plan's results", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And no plan named \"BDD plan\" is saved", () => noSavedPlan(page, "BDD plan"));
    await session.step(9, "And user opens PMX Studio", () => openApp(page));
    await session.step(10, "Then the study should have finished running", () => finishedRunning(page));
    await run.scenario("Design tab", async () => {
      await session.step(13, "Then \"Experimental groups\" section should contain text \"Treatment 1\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Treatment 1"));
      await session.step(14, "And \"Predicted tumor growth\" section should be visible", () => shouldBe(page, el("\"Predicted tumor growth\" section"), "visible"));
      await session.step(15, "And visible volume axis should have the value \"Log\"", () => shouldHaveValue(page, el("visible volume axis"), "Log"));
      await session.step(16, "When user clicks on linear axis choice", () => clickOn(page, el("linear axis choice")));
      await session.step(17, "Then visible volume axis should have the value \"Linear\"", () => shouldHaveValue(page, el("visible volume axis"), "Linear"));
      await session.step(18, "When user clicks on log axis choice", () => clickOn(page, el("log axis choice")));
      await session.step(19, "Then visible volume axis should have the value \"Log\"", () => shouldHaveValue(page, el("visible volume axis"), "Log"));
    });
    await run.scenario("Calendar tab, measurement days", async () => {
      await session.step(22, "When user clicks on \"Calendar\" tab", () => clickOn(page, el("\"Calendar\" tab")));
      await session.step(23, "Then \"Study calendar\" section should be visible", () => shouldBe(page, el("\"Study calendar\" section"), "visible"));
      await session.step(24, "And calendar summary should contain text \"9 measurement days\"", () => shouldContainText(page, el("calendar summary"), "9 measurement days"));
      await session.step(25, "When user clicks on tuesday header", () => clickOn(page, el("tuesday header")));
      await session.step(26, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(27, "And calendar summary should contain text \"12 measurement days\"", () => shouldContainText(page, el("calendar summary"), "12 measurement days"));
      await session.step(28, "When user clicks on tuesday header", () => clickOn(page, el("tuesday header")));
      await session.step(29, "Then calendar summary should contain text \"9 measurement days\"", () => shouldContainText(page, el("calendar summary"), "9 measurement days"));
    });
    await run.scenario("Calendar tab, randomization day", async () => {
      await session.step(32, "When user selects \"Wednesday\" in randomization day picker", () => selectIn(page, "Wednesday", el("randomization day picker")));
      await session.step(33, "Then calendar summary should contain text \"Randomize on a\"", () => shouldContainText(page, el("calendar summary"), "Randomize on a"));
      await session.step(34, "And randomization day picker should have the value \"Wednesday\"", () => shouldHaveValue(page, el("randomization day picker"), "Wednesday"));
      await session.step(35, "When user selects \"Monday\" in randomization day picker", () => selectIn(page, "Monday", el("randomization day picker")));
      await session.step(36, "Then randomization day picker should have the value \"Monday\"", () => shouldHaveValue(page, el("randomization day picker"), "Monday"));
      await session.step(37, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(38, "Then \"Study protocol\" section should contain text \"days 0, 2, 4, 7\"", () => shouldContainText(page, el("\"Study protocol\" section"), "days 0, 2, 4, 7"));
      await session.step(39, "And randomization day picker should have the value \"Monday\"", () => shouldHaveValue(page, el("randomization day picker"), "Monday"));
    });
    await run.scenario("Calendar tab, two changes without waiting are both planned", async () => {
      await session.step(42, "When user clicks on \"Calendar\" tab", () => clickOn(page, el("\"Calendar\" tab")));
      await session.step(43, "And user clicks on tuesday header", () => clickOn(page, el("tuesday header")));
      await session.step(44, "And user clicks on thursday header", () => clickOn(page, el("thursday header")));
      await session.step(45, "And user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(46, "Then \"Study protocol\" section should contain text \"Mon, Tue, Wed, Thu, Fri\"", () => shouldContainText(page, el("\"Study protocol\" section"), "Mon, Tue, Wed, Thu, Fri"));
      await session.step(47, "When user clicks on \"Calendar\" tab", () => clickOn(page, el("\"Calendar\" tab")));
      await session.step(48, "And user clicks on tuesday header", () => clickOn(page, el("tuesday header")));
      await session.step(49, "And user clicks on thursday header", () => clickOn(page, el("thursday header")));
      await session.step(50, "And user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(51, "Then \"Study protocol\" section should contain text \"on Mon, Wed, Fri\"", () => shouldContainText(page, el("\"Study protocol\" section"), "on Mon, Wed, Fri"));
    });
    await run.scenario("A form edited during a replan stays out of date", async () => {
      await session.step(54, "When user clicks on \"Calendar\" tab", () => clickOn(page, el("\"Calendar\" tab")));
      await session.step(55, "And user selects \"Tuesday\" in randomization day picker", () => selectIn(page, "Tuesday", el("randomization day picker")));
      await session.step(56, "And user enters \"1400\" into \"Endpoint volume\" input", () => enterInto(page, "1400", el("\"Endpoint volume\" input")));
      await session.step(57, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(58, "And randomization day picker should have the value \"Tuesday\"", () => shouldHaveValue(page, el("randomization day picker"), "Tuesday"));
      await session.step(59, "And stale note should be visible", () => shouldBe(page, el("stale note"), "visible"));
      await session.step(60, "When user enters \"1500\" into \"Endpoint volume\" input", () => enterInto(page, "1500", el("\"Endpoint volume\" input")));
      await session.step(61, "And user selects \"Monday\" in randomization day picker", () => selectIn(page, "Monday", el("randomization day picker")));
      await session.step(62, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(63, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(64, "And stale note should be hidden", () => shouldBe(page, el("stale note"), "hidden"));
    });
    await run.scenario("Trade-offs tab", async () => {
      await session.step(67, "When user clicks on \"Trade-offs\" tab", () => clickOn(page, el("\"Trade-offs\" tab")));
      await session.step(68, "Then \"Power by group size\" section should be visible", () => shouldBe(page, el("\"Power by group size\" section"), "visible"));
      await session.step(69, "And \"Power by end of study\" section should contain text \"Ending after day\"", () => shouldContainText(page, el("\"Power by end of study\" section"), "Ending after day"));
      await session.step(70, "And \"Measurement schedule\" section should contain text \"this plan\"", () => shouldContainText(page, el("\"Measurement schedule\" section"), "this plan"));
      await session.step(71, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(72, "And user enters \"1400\" into \"Endpoint volume\" input", () => enterInto(page, "1400", el("\"Endpoint volume\" input")));
      await session.step(73, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(74, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(75, "When user enters \"1500\" into \"Endpoint volume\" input", () => enterInto(page, "1500", el("\"Endpoint volume\" input")));
      await session.step(76, "And user clicks on \"Trade-offs\" tab", () => clickOn(page, el("\"Trade-offs\" tab")));
      await session.step(77, "Then \"Measurement schedule\" section should contain text \"Run to compare\"", () => shouldContainText(page, el("\"Measurement schedule\" section"), "Run to compare"));
      await session.step(78, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(79, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(80, "And \"Measurement schedule\" section should contain text \"this plan\"", () => shouldContainText(page, el("\"Measurement schedule\" section"), "this plan"));
    });
    await run.scenario("Protocol tab", async () => {
      await session.step(83, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(84, "Then \"Study protocol\" section should contain text \"MC38\"", () => shouldContainText(page, el("\"Study protocol\" section"), "MC38"));
      await session.step(85, "And \"Study protocol\" section should contain text \"to implant\"", () => shouldContainText(page, el("\"Study protocol\" section"), "to implant"));
      await session.step(86, "When user clicks on \"Copy\" button", () => clickOn(page, el("\"Copy\" button")));
      await session.step(87, "Then an info balloon containing \"Study protocol copied\" should have been shown", () => infoBalloonText(page, "Study protocol copied"));
      await session.step(88, "And the clipboard should contain text \"Study protocol (Plan MC38\"", () => clipboardContains(page, "Study protocol (Plan MC38"));
    });
    await run.scenario("Saving a plan, and replacing it on a second press", async () => {
      await session.step(91, "When user clicks on \"Save plan\" button", () => clickOn(page, el("\"Save plan\" button")));
      await session.step(92, "Then \"Save plan\" dialog should be visible", () => shouldBe(page, el("\"Save plan\" dialog"), "visible"));
      await session.step(93, "When user enters \"BDD plan\" into \"Name\" input in \"Save plan\" dialog", () => enterInto(page, "BDD plan", el("\"Name\" input in \"Save plan\" dialog")));
      await session.step(94, "And user clicks on \"Save\" button in \"Save plan\" dialog", () => clickOn(page, el("\"Save\" button in \"Save plan\" dialog")));
      await session.step(95, "Then an info balloon containing \"Plan saved\" should have been shown", () => infoBalloonText(page, "Plan saved"));
      await session.step(96, "And a plan named \"BDD plan\" should be saved", () => planIsSaved(page, "BDD plan"));
      await session.step(97, "And \"Save plan\" dialog should be absent", () => shouldBe(page, el("\"Save plan\" dialog"), "absent"));
      await session.step(98, "When user clicks on \"Save plan\" button", () => clickOn(page, el("\"Save plan\" button")));
      await session.step(99, "And user enters \"BDD plan\" into \"Name\" input in \"Save plan\" dialog", () => enterInto(page, "BDD plan", el("\"Name\" input in \"Save plan\" dialog")));
      await session.step(100, "And user clicks on \"Save\" button in \"Save plan\" dialog", () => clickOn(page, el("\"Save\" button in \"Save plan\" dialog")));
      await session.step(101, "Then \"Save plan\" dialog should contain text \"already saved\"", () => shouldContainText(page, el("\"Save plan\" dialog"), "already saved"));
      await session.step(102, "When user clicks on \"Save\" button in \"Save plan\" dialog", () => clickOn(page, el("\"Save\" button in \"Save plan\" dialog")));
      await session.step(103, "Then \"Save plan\" dialog should be absent", () => shouldBe(page, el("\"Save plan\" dialog"), "absent"));
    });
    await run.scenario("Downloading the plan", async () => {
      await session.step(106, "Given user watches downloads", () => watchDownloads(page));
      await session.step(107, "When user downloads a file through \"Download\" button", () => downloadThrough(page, el("\"Download\" button")));
      await session.step(108, "Then a file matching \"\\.pmxplan$\" should have been downloaded", () => fileMatchingDownloaded(page, "\\.pmxplan$"));
      await session.step(109, "And the downloaded file should contain \"pmx-plan\"", () => downloadedContains(page, "pmx-plan"));
    });
    await run.scenario("A plan out of date is not saved", async () => {
      await session.step(112, "When user enters \"1200\" into \"Endpoint volume\" input", () => enterInto(page, "1200", el("\"Endpoint volume\" input")));
      await session.step(113, "And user clicks on \"Save plan\" button", () => clickOn(page, el("\"Save plan\" button")));
      await session.step(114, "Then a warning balloon containing \"run to update the plan first\" should have been shown", () => warningBalloonText(page, "run to update the plan first"));
      await session.step(115, "And \"Save plan\" dialog should be absent", () => shouldBe(page, el("\"Save plan\" dialog"), "absent"));
      await session.step(116, "When user enters \"1500\" into \"Endpoint volume\" input", () => enterInto(page, "1500", el("\"Endpoint volume\" input")));
      await session.step(117, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(118, "Then the study should have finished running", () => finishedRunning(page));
    });
    await run.scenario("Moving between the stages", async () => {
      await session.step(121, "When user clicks on \"3 Analyze\" button", () => clickOn(page, el("\"3 Analyze\" button")));
      await session.step(122, "Then analyze start should be visible", () => shouldBe(page, el("analyze start"), "visible"));
      await session.step(123, "When user clicks on \"1 Plan\" button", () => clickOn(page, el("\"1 Plan\" button")));
      await session.step(124, "Then \"Study protocol\" section should be visible", () => shouldBe(page, el("\"Study protocol\" section"), "visible"));
    });
    run.finish();
  });
});
