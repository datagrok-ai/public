/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/forms/func-calls.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.func-form]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {clearField, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveError} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("FuncCalls", () => {
  const session = feature(test, "features/demo/forms/func-calls.feature", import.meta.url);
  test("A stage picked on either side writes the call and shows on the other", {tag: ["@demo", "@realizes:u2.func-form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the \"FuncCalls\" demo page", () => openDemoPage(page, "FuncCalls"));
    enter(page, "U2 Demo");
    await session.step(12, "When user selects \"Phase II\" in Stage input in \"u2\" section", () => selectIn(page, "Phase II", el("Stage input in \"u2\" section")));
    await session.step(13, "Then Stage input in \"Platform\" section should have value \"Phase II\"", () => shouldHaveValue(page, el("Stage input in \"Platform\" section"), "Phase II"));
    await session.step(14, "And value of first inputs readout should contain text \"stage: 'Phase II'\"", () => shouldContainText(page, el("value of first inputs readout"), "stage: 'Phase II'"));
    await session.step(15, "And value of first sync readout should have text \"u2 → Dart 1 · Dart → u2 0\"", () => shouldHaveText(page, el("value of first sync readout"), "u2 → Dart 1 · Dart → u2 0"));
    await session.step(16, "When user selects \"Discovery\" in Stage input in \"Platform\" section", () => selectIn(page, "Discovery", el("Stage input in \"Platform\" section")));
    await session.step(17, "Then Stage input in \"u2\" section should have value \"Discovery\"", () => shouldHaveValue(page, el("Stage input in \"u2\" section"), "Discovery"));
    await session.step(18, "And value of first inputs readout should contain text \"stage: 'Discovery'\"", () => shouldContainText(page, el("value of first inputs readout"), "stage: 'Discovery'"));
    await session.step(19, "And value of first sync readout should have text \"u2 → Dart 1 · Dart → u2 1\"", () => shouldHaveText(page, el("value of first sync readout"), "u2 → Dart 1 · Dart → u2 1"));
  });
  test("Emptying the required parameter makes both forms invalid", {tag: ["@demo", "@realizes:u2.func-form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the \"FuncCalls\" demo page", () => openDemoPage(page, "FuncCalls"));
    enter(page, "U2 Demo");
    await session.step(22, "Then value of isValid readout should have text \"Dart true · u2 false\"", () => shouldHaveText(page, el("value of isValid readout"), "Dart true · u2 false"));
    await session.step(23, "When user enters \"B-17\" into \"Batch Id\" input in \"u2\" section", () => enterInto(page, "B-17", el("\"Batch Id\" input in \"u2\" section")));
    await session.step(24, "Then \"Batch Id\" input in \"Platform\" section should have value \"B-17\"", () => shouldHaveValue(page, el("\"Batch Id\" input in \"Platform\" section"), "B-17"));
    await session.step(25, "And \"Batch Id\" input in \"Platform\" section should be valid", () => shouldBe(page, el("\"Batch Id\" input in \"Platform\" section"), "valid"));
    await session.step(26, "And value of isValid readout should have text \"Dart true · u2 true\"", () => shouldHaveText(page, el("value of isValid readout"), "Dart true · u2 true"));
    await session.step(27, "When user clears \"Batch Id\" input in \"u2\" section", () => clearField(page, el("\"Batch Id\" input in \"u2\" section")));
    await session.step(28, "Then \"Batch Id\" input in \"Platform\" section should be invalid", () => shouldBe(page, el("\"Batch Id\" input in \"Platform\" section"), "invalid"));
    await session.step(29, "And \"Batch Id\" input in \"u2\" section should have the error \"Value can't be empty\"", () => shouldHaveError(page, el("\"Batch Id\" input in \"u2\" section"), "Value can't be empty"));
    await session.step(30, "And value of isValid readout should have text \"Dart false · u2 false\"", () => shouldHaveText(page, el("value of isValid readout"), "Dart false · u2 false"));
    await session.step(31, "When user enters \"B-18\" into \"Batch Id\" input in \"Platform\" section", () => enterInto(page, "B-18", el("\"Batch Id\" input in \"Platform\" section")));
    await session.step(32, "Then \"Batch Id\" input in \"u2\" section should have value \"B-18\"", () => shouldHaveValue(page, el("\"Batch Id\" input in \"u2\" section"), "B-18"));
    await session.step(33, "And \"Batch Id\" input in \"u2\" section should be valid", () => shouldBe(page, el("\"Batch Id\" input in \"u2\" section"), "valid"));
    await session.step(34, "And value of isValid readout should have text \"Dart true · u2 true\"", () => shouldHaveText(page, el("value of isValid readout"), "Dart true · u2 true"));
  });
});
