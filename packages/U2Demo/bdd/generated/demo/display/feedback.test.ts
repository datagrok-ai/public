/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/display/feedback.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.progress, u2.notify, u2.tour]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, close, hoverOver, pressKey, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Progress, notifications and the guided tour", () => {
  const session = feature(test, "features/demo/display/feedback.feature", import.meta.url);
  test("Progress runs to completion and resets", {tag: ["@demo", "@realizes:u2.progress", "@realizes:u2.notify", "@realizes:u2.tour"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(5, "Given user opens the \"Feedback\" demo page", () => openDemoPage(page, "Feedback"));
    enter(page, "U2 Demo");
    await session.step(8, "Then \"Indexing compounds\" progress bar should contain text \"0%\"", () => shouldContainText(page, el("\"Indexing compounds\" progress bar"), "0%"));
    await session.step(9, "When user clicks on Run button", () => clickOn(page, el("Run button")));
    await session.step(10, "Then \"Indexing compounds\" progress bar should contain text \"100%\"", () => shouldContainText(page, el("\"Indexing compounds\" progress bar"), "100%"));
    await session.step(11, "And value of progress readout should have text \"1.00\"", () => shouldHaveText(page, el("value of progress readout"), "1.00"));
    await session.step(12, "When user clicks on Reset button", () => clickOn(page, el("Reset button")));
    await session.step(13, "Then value of progress readout should have text \"0.00\"", () => shouldHaveText(page, el("value of progress readout"), "0.00"));
  });
  test("Notifications stack and close", {tag: ["@demo", "@realizes:u2.progress", "@realizes:u2.notify", "@realizes:u2.tour"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(5, "Given user opens the \"Feedback\" demo page", () => openDemoPage(page, "Feedback"));
    enter(page, "U2 Demo");
    await session.step(16, "When user clicks on Info button", () => clickOn(page, el("Info button")));
    await session.step(17, "Then notification should contain text \"Indexed 1,204 compounds.\"", () => shouldContainText(page, el("notification"), "Indexed 1,204 compounds."));
    await session.step(18, "When user clicks on Warning button", () => clickOn(page, el("Warning button")));
    await session.step(19, "And user clicks on Error button", () => clickOn(page, el("Error button")));
    await session.step(20, "Then \"Two structures could not be parsed.\" notification should be visible", () => shouldBe(page, el("\"Two structures could not be parsed.\" notification"), "visible"));
    await session.step(21, "And \"Connection refused.\" notification should be visible", () => shouldBe(page, el("\"Connection refused.\" notification"), "visible"));
    await session.step(22, "When user hovers over \"Connection refused.\" notification", () => hoverOver(page, el("\"Connection refused.\" notification")));
    await session.step(23, "And user closes \"Connection refused.\" notification", () => close(page, el("\"Connection refused.\" notification")));
    await session.step(24, "Then \"Connection refused.\" notification should be hidden", () => shouldBe(page, el("\"Connection refused.\" notification"), "hidden"));
    await session.step(25, "When user clicks on \"Close all\" button", () => clickOn(page, el("\"Close all\" button")));
    await session.step(26, "Then \"Two structures could not be parsed.\" notification should be hidden", () => shouldBe(page, el("\"Two structures could not be parsed.\" notification"), "hidden"));
  });
  test("The tour runs to its end and skips a missing target", {tag: ["@demo", "@realizes:u2.progress", "@realizes:u2.notify", "@realizes:u2.tour"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(5, "Given user opens the \"Feedback\" demo page", () => openDemoPage(page, "Feedback"));
    enter(page, "U2 Demo");
    await session.step(29, "When user clicks on \"Start tour\" button", () => clickOn(page, el("\"Start tour\" button")));
    await session.step(30, "Then tour should be visible", () => shouldBe(page, el("tour"), "visible"));
    await session.step(31, "And \"1 / 3\" text should be visible", () => shouldBe(page, el("\"1 / 3\" text"), "visible"));
    await session.step(32, "When user clicks on NEXT button", () => clickOn(page, el("NEXT button")));
    await session.step(33, "And user clicks on NEXT button", () => clickOn(page, el("NEXT button")));
    await session.step(34, "And user clicks on DONE button", () => clickOn(page, el("DONE button")));
    await session.step(35, "Then tour should be hidden", () => shouldBe(page, el("tour"), "hidden"));
    await session.step(36, "And value of tour readout should have text \"done\"", () => shouldHaveText(page, el("value of tour readout"), "done"));
    await session.step(37, "When user clicks on \"Start tour\" button", () => clickOn(page, el("\"Start tour\" button")));
    await session.step(38, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(39, "Then value of tour readout should have text \"skipped\"", () => shouldHaveText(page, el("value of tour readout"), "skipped"));
  });
});
