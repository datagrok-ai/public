/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/inputs/chained-choices.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.inputs]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {selectIn, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldBeDisabledBecause, shouldListInOrder} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Chained choices", () => {
  const session = feature(test, "features/demo/inputs/chained-choices.feature", import.meta.url);
  test("Each pick loads the next choice's options", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Chained choices\" demo page", () => openDemoPage(page, "Chained choices"));
    enter(page, "U2 Demo");
    await session.step(11, "Then Letter input should be disabled because \"Pick Name first\"", () => shouldBeDisabledBecause(page, el("Letter input"), "Pick Name first"));
    await session.step(12, "And Shape input should be disabled because \"Pick Name first\"", () => shouldBeDisabledBecause(page, el("Shape input"), "Pick Name first"));
    await session.step(13, "When user selects \"Ada\" in Name input", () => selectIn(page, "Ada", el("Name input")));
    await session.step(14, "And user selects \"Ada a\" in Letter input", () => selectIn(page, "Ada a", el("Letter input")));
    await session.step(15, "And user selects \"Ada a 1\" in Digit input", () => selectIn(page, "Ada a 1", el("Digit input")));
    await session.step(16, "Then value of picked readout should have text \"Ada → Ada a → Ada a 1 → – → –\"", () => shouldHaveText(page, el("value of picked readout"), "Ada → Ada a → Ada a 1 → – → –"));
    await session.step(17, "And Sign input should list the items \"Ada a 1 +, Ada a 1 =\" in this order", () => shouldListInOrder(page, el("Sign input"), "Ada a 1 +, Ada a 1 ="));
    await session.step(18, "And Shape input should be disabled because \"Pick Sign first\"", () => shouldBeDisabledBecause(page, el("Shape input"), "Pick Sign first"));
  });
  test("Changing a pick clears what was built on it and reloads the lists below", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Chained choices\" demo page", () => openDemoPage(page, "Chained choices"));
    enter(page, "U2 Demo");
    await session.step(21, "When user selects \"Ada\" in Name input", () => selectIn(page, "Ada", el("Name input")));
    await session.step(22, "And user selects \"Ada b\" in Letter input", () => selectIn(page, "Ada b", el("Letter input")));
    await session.step(23, "And user selects \"Ada b 3\" in Digit input", () => selectIn(page, "Ada b 3", el("Digit input")));
    await session.step(24, "And user selects \"Cleo\" in Name input", () => selectIn(page, "Cleo", el("Name input")));
    await session.step(25, "Then value of picked readout should have text \"Cleo → – → – → – → –\"", () => shouldHaveText(page, el("value of picked readout"), "Cleo → – → – → – → –"));
    await session.step(26, "And Letter input should list the items \"Cleo a, Cleo e\" in this order", () => shouldListInOrder(page, el("Letter input"), "Cleo a, Cleo e"));
    await session.step(27, "And Digit input should be disabled because \"Pick Letter first\"", () => shouldBeDisabledBecause(page, el("Digit input"), "Pick Letter first"));
  });
});
