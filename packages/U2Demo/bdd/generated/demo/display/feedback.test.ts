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
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, close, pressKey, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Progress, notifications and the guided tour", () => {
  const session = feature(test);
  test("Progress runs to completion and resets", {tag: ["@demo", "@realizes:u2.progress", "@realizes:u2.notify", "@realizes:u2.tour"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Feedback\" demo page", () => openDemoPage(page, "Feedback"));
    enter(page, "U2 Demo");
    await test.step("Then \"Indexing compounds\" progress bar should contain text \"0%\"", () => shouldContainText(page, el("\"Indexing compounds\" progress bar"), "0%"));
    await test.step("When user clicks on Run button", () => clickOn(page, el("Run button")));
    await test.step("Then \"Indexing compounds\" progress bar should contain text \"100%\"", () => shouldContainText(page, el("\"Indexing compounds\" progress bar"), "100%"));
    await test.step("And value of progress readout should have text \"1.00\"", () => shouldHaveText(page, el("value of progress readout"), "1.00"));
    await test.step("When user clicks on Reset button", () => clickOn(page, el("Reset button")));
    await test.step("Then value of progress readout should have text \"0.00\"", () => shouldHaveText(page, el("value of progress readout"), "0.00"));
  });
  test("Notifications stack and close", {tag: ["@demo", "@realizes:u2.progress", "@realizes:u2.notify", "@realizes:u2.tour"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Feedback\" demo page", () => openDemoPage(page, "Feedback"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on Info button", () => clickOn(page, el("Info button")));
    await test.step("Then notification should contain text \"Indexed 1,204 compounds.\"", () => shouldContainText(page, el("notification"), "Indexed 1,204 compounds."));
    await test.step("When user clicks on Warning button", () => clickOn(page, el("Warning button")));
    await test.step("And user clicks on Error button", () => clickOn(page, el("Error button")));
    await test.step("Then \"Two structures could not be parsed.\" notification should be visible", () => shouldBe(page, el("\"Two structures could not be parsed.\" notification"), "visible"));
    await test.step("And \"Connection refused.\" notification should be visible", () => shouldBe(page, el("\"Connection refused.\" notification"), "visible"));
    await test.step("When user closes \"Connection refused.\" notification", () => close(page, el("\"Connection refused.\" notification")));
    await test.step("Then \"Connection refused.\" notification should be hidden", () => shouldBe(page, el("\"Connection refused.\" notification"), "hidden"));
    await test.step("When user clicks on \"Close all\" button", () => clickOn(page, el("\"Close all\" button")));
    await test.step("Then \"Two structures could not be parsed.\" notification should be hidden", () => shouldBe(page, el("\"Two structures could not be parsed.\" notification"), "hidden"));
  });
  test("The tour runs to its end and skips a missing target", {tag: ["@demo", "@realizes:u2.progress", "@realizes:u2.notify", "@realizes:u2.tour"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Feedback\" demo page", () => openDemoPage(page, "Feedback"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"Start tour\" button", () => clickOn(page, el("\"Start tour\" button")));
    await test.step("Then tour should be visible", () => shouldBe(page, el("tour"), "visible"));
    await test.step("And \"1 / 4\" text should be visible", () => shouldBe(page, el("\"1 / 4\" text"), "visible"));
    await test.step("When user clicks on NEXT button", () => clickOn(page, el("NEXT button")));
    await test.step("And user clicks on NEXT button", () => clickOn(page, el("NEXT button")));
    await test.step("And user clicks on NEXT button", () => clickOn(page, el("NEXT button")));
    await test.step("Then tour should be hidden", () => shouldBe(page, el("tour"), "hidden"));
    await test.step("And value of tour readout should have text \"done\"", () => shouldHaveText(page, el("value of tour readout"), "done"));
    await test.step("When user clicks on \"Start tour\" button", () => clickOn(page, el("\"Start tour\" button")));
    await test.step("And user presses Escape", () => pressKey(page, "Escape"));
    await test.step("Then value of tour readout should have text \"skipped\"", () => shouldHaveText(page, el("value of tour readout"), "skipped"));
  });
});
