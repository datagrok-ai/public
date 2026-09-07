/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/bridge.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dg.dart-input]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, shouldBe, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Dart bridge", () => {
  const session = feature(test, "features/demo/platform/bridge.feature", import.meta.url);
  test("A u2 input inside a platform dialog", {tag: ["@demo", "@realizes:u2.dg.dart-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Bridge\" demo page", () => openDemoPage(page, "Bridge"));
    enter(page, "U2 Demo");
    await session.step(8, "When user clicks on \"Open a DG.Dialog with a u2 input\" button", () => clickOn(page, el("\"Open a DG.Dialog with a u2 input\" button")));
    await session.step(9, "Then dialog should be visible", () => shouldBe(page, el("dialog"), "visible"));
    await session.step(10, "And title of dialog should have text \"u2 input, platform dialog\"", () => shouldHaveText(page, el("title of dialog"), "u2 input, platform dialog"));
    await session.step(11, "When user types \"Caffeine\" into Compound input in dialog", () => typeInto(page, "Caffeine", el("Compound input in dialog")));
    await session.step(12, "Then value of \"bridged value\" readout should have text \"Caffeine\"", () => shouldHaveText(page, el("value of \"bridged value\" readout"), "Caffeine"));
    await session.step(13, "When user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await session.step(14, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
    await session.step(15, "And notification should contain text \"Compound: Caffeine\"", () => shouldContainText(page, el("notification"), "Compound: Caffeine"));
  });
  test("The leak detector answers", {tag: ["@demo", "@realizes:u2.dg.dart-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user opens the \"Bridge\" demo page", () => openDemoPage(page, "Bridge"));
    enter(page, "U2 Demo");
    await session.step(19, "When user clicks on \"leakReport()\" button", () => clickOn(page, el("\"leakReport()\" button")));
    await session.step(20, "Then demo page should contain text \"liveScopes\"", () => shouldContainText(page, el("demo page"), "liveScopes"));
  });
});
