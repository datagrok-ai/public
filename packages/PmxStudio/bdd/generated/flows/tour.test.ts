/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/flows/tour.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {finishedRunning} from '../../bindings/pmx.js';
import {openApp} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragBy, dragSliderTo, dragTo, scrollIntoView, selectIn, shouldBe, shouldContainText, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A tour of PMX Studio", () => {
  const session = feature(test, "features/flows/tour.feature", import.meta.url);
  test("A tour of PMX Studio", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await run.scenario("From the Overview to the analysis of an example study", async () => {
      await session.step(9, "Given user is logged in", () => loggedIn(page));
      await session.step(10, "And user opens PMX Studio", () => openApp(page));
      await session.step(12, "When user clicks on plan card", () => clickOn(page, el("plan card")), undefined, "Start a new study in Plan");
      await session.step(14, "And user clicks on \"Add treatment arm\" button", () => clickOn(page, el("\"Add treatment arm\" button")), undefined, "A new plan starts with no arm: add one");
      await session.step(15, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(16, "And \"Experimental groups\" section should be visible", () => shouldBe(page, el("\"Experimental groups\" section"), "visible"));
      await session.step(18, "When user clicks on overview step", () => clickOn(page, el("overview step")), undefined, "Back to the Overview: open the example study, a table of measurements");
      await session.step(19, "And user clicks on open a study strip", () => clickOn(page, el("open a study strip")));
      await session.step(20, "And user uploads \"../files/example/mc38-ld4172-anti-pd1.csv\" through study drop zone", () => uploadThrough(page, "../files/example/mc38-ld4172-anti-pd1.csv", el("study drop zone")));
      await session.step(22, "And user clicks on \"Open in Analyze\" button in \"Open a study\" dialog", () => clickOn(page, el("\"Open in Analyze\" button in \"Open a study\" dialog")), undefined, "It has finished: analyze it");
      await session.step(23, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(24, "And \"Observed tumor growth\" section should be visible", () => shouldBe(page, el("\"Observed tumor growth\" section"), "visible"));
      await session.step(26, "When user picks \"View > Layout…\" from the top menu", () => pickFromTopMenu(page, "View > Layout…"), undefined, "Arrange the panels: View > Layout");
      await session.step(27, "Then \"Layout of Analyze · Results\" dialog should be visible", () => shouldBe(page, el("\"Layout of Analyze · Results\" dialog"), "visible"));
      await session.step(29, "When user drags growth rates corner by 120 and 0 pixels", () => dragBy(page, el("growth rates corner"), 120, 0), undefined, "Widen a panel from its corner");
      await session.step(31, "And user drags survival tile to growth rates tile", () => dragTo(page, el("survival tile"), el("growth rates tile")), undefined, "Move a panel: the others make room");
      await session.step(33, "And user clicks on \"OK\" button in \"Layout of Analyze · Results\" dialog", () => clickOn(page, el("\"OK\" button in \"Layout of Analyze · Results\" dialog")), undefined, "Keep the layout");
      await session.step(35, "When user clicks on each animal tab", () => clickOn(page, el("each animal tab")), undefined, "Every animal, compared by what you choose");
      await session.step(36, "And user scrolls waterfall chart into view", () => scrollIntoView(page, el("waterfall chart")));
      await session.step(37, "And user selects \"Best response\" in \"Compare by\" input", () => selectIn(page, "Best response", el("\"Compare by\" input")));
      await session.step(39, "When user scrolls first animal row into view", () => scrollIntoView(page, el("first animal row")), undefined, "Pick an animal: every chart marks it");
      await session.step(40, "And user clicks on first animal row", () => clickOn(page, el("first animal row")));
      await session.step(42, "When user clicks on results tab", () => clickOn(page, el("results tab")), undefined, "Follow the study past its last measurement");
      await session.step(43, "And user scrolls km chart into view", () => scrollIntoView(page, el("km chart")));
      await session.step(44, "And user drags the slider of km slider to 24", () => dragSliderTo(page, el("km slider"), 24));
      await session.step(45, "Then \"Time to endpoint\" section should contain text \"projected\"", () => shouldContainText(page, el("\"Time to endpoint\" section"), "projected"));
    });
    run.finish();
  });
});
