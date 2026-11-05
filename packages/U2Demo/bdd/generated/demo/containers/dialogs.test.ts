/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/containers/dialogs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dialog]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {tallPx} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {clickOn, pressKey, shouldBe, shouldContainText, shouldHaveText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Dialogs and confirm", () => {
  const session = feature(test, "features/demo/containers/dialogs.feature", import.meta.url);
  test("Delete asks first; CANCEL deletes nothing, DELETE deletes once its work is done", {tag: ["@demo", "@realizes:u2.dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Dialogs\" demo page", () => openDemoPage(page, "Dialogs"));
    enter(page, "U2 Demo");
    await session.step(11, "When user clicks on \"Delete campaign...\" button", () => clickOn(page, el("\"Delete campaign...\" button")));
    await session.step(12, "Then \"Delete campaign\" dialog should be visible", () => shouldBe(page, el("\"Delete campaign\" dialog"), "visible"));
    await session.step(13, "And \"Delete campaign\" dialog should contain text \"This cannot be undone.\"", () => shouldContainText(page, el("\"Delete campaign\" dialog"), "This cannot be undone."));
    await session.step(14, "When user clicks on CANCEL button in \"Delete campaign\" dialog", () => clickOn(page, el("CANCEL button in \"Delete campaign\" dialog")));
    await session.step(15, "Then \"Delete campaign\" dialog should be hidden", () => shouldBe(page, el("\"Delete campaign\" dialog"), "hidden"));
    await session.step(16, "And status bar should contain text \"Delete cancelled\"", () => shouldContainText(page, el("status bar"), "Delete cancelled"));
    await session.step(17, "And value of deleted readout should have text \"0\"", () => shouldHaveText(page, el("value of deleted readout"), "0"));
    await session.step(18, "When user clicks on \"Delete campaign...\" button", () => clickOn(page, el("\"Delete campaign...\" button")));
    await session.step(19, "And user clicks on Delete button in \"Delete campaign\" dialog", () => clickOn(page, el("Delete button in \"Delete campaign\" dialog")));
    await session.step(20, "Then \"Delete campaign\" dialog should be hidden", () => shouldBe(page, el("\"Delete campaign\" dialog"), "hidden"));
    await session.step(21, "And status bar should contain text \"Campaign KHL-3 deleted\"", () => shouldContainText(page, el("status bar"), "Campaign KHL-3 deleted"));
    await session.step(22, "And value of deleted readout should have text \"1\"", () => shouldHaveText(page, el("value of deleted readout"), "1"));
  });
  test("The first focus is CANCEL, never the action, so Enter deletes nothing", {tag: ["@demo", "@realizes:u2.dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Dialogs\" demo page", () => openDemoPage(page, "Dialogs"));
    enter(page, "U2 Demo");
    await session.step(25, "When user clicks on \"Delete campaign...\" button", () => clickOn(page, el("\"Delete campaign...\" button")));
    await session.step(26, "Then CANCEL button in \"Delete campaign\" dialog should be focused", () => shouldBe(page, el("CANCEL button in \"Delete campaign\" dialog"), "focused"));
    await session.step(27, "When user presses Enter", () => pressKey(page, "Enter"));
    await session.step(28, "Then \"Delete campaign\" dialog should be hidden", () => shouldBe(page, el("\"Delete campaign\" dialog"), "hidden"));
    await session.step(29, "And status bar should contain text \"Delete cancelled\"", () => shouldContainText(page, el("status bar"), "Delete cancelled"));
    await session.step(30, "And value of deleted readout should have text \"0\"", () => shouldHaveText(page, el("value of deleted readout"), "0"));
  });
  test("A fixed height: the footer stays at the bottom", {tag: ["@demo", "@realizes:u2.dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Dialogs\" demo page", () => openDemoPage(page, "Dialogs"));
    enter(page, "U2 Demo");
    await session.step(33, "When user clicks on \"Sized dialog\" button", () => clickOn(page, el("\"Sized dialog\" button")));
    await session.step(34, "Then \"Progress tracker\" dialog should be 420 px tall", () => tallPx(page, el("\"Progress tracker\" dialog"), 420));
    await session.step(35, "And body of \"Progress tracker\" dialog should be 339 px tall", () => tallPx(page, el("body of \"Progress tracker\" dialog"), 339));
    await session.step(36, "When user clicks on CANCEL button in \"Progress tracker\" dialog", () => clickOn(page, el("CANCEL button in \"Progress tracker\" dialog")));
    await session.step(37, "Then \"Progress tracker\" dialog should be hidden", () => shouldBe(page, el("\"Progress tracker\" dialog"), "hidden"));
  });
  test("The first focus goes to the first field that takes edits", {tag: ["@demo", "@realizes:u2.dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Dialogs\" demo page", () => openDemoPage(page, "Dialogs"));
    enter(page, "U2 Demo");
    await session.step(40, "When user clicks on \"Submit...\" button", () => clickOn(page, el("\"Submit...\" button")));
    await session.step(41, "Then Status input in Submit dialog should be focused", () => shouldBe(page, el("Status input in Submit dialog"), "focused"));
    await session.step(42, "And Name input in Submit dialog should be read-only", () => shouldBe(page, el("Name input in Submit dialog"), "read-only"));
    await session.step(43, "When user clicks on CANCEL button in Submit dialog", () => clickOn(page, el("CANCEL button in Submit dialog")));
    await session.step(44, "Then Submit dialog should be hidden", () => shouldBe(page, el("Submit dialog"), "hidden"));
  });
  test("Required fields say so", {tag: ["@demo", "@realizes:u2.dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Dialogs\" demo page", () => openDemoPage(page, "Dialogs"));
    enter(page, "U2 Demo");
    await session.step(47, "Then \"template-name\" text input should be required", () => shouldBe(page, el("\"template-name\" text input"), "required"));
    await session.step(48, "And \"template-key\" text input should be required", () => shouldBe(page, el("\"template-key\" text input"), "required"));
    await session.step(49, "And \"template-description\" text input should not be required", () => shouldNotBe(page, el("\"template-description\" text input"), "required"));
  });
});
