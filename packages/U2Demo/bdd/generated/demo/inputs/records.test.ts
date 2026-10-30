/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/inputs/records.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.records-input, u2.checklist-input]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, pressKeyIn, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {itemShouldBe, shouldHaveItemsInTotal, shouldListInOrder} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Records and an ordered checklist", () => {
  const session = feature(test, "features/demo/inputs/records.feature", import.meta.url);
  test("A new row is named in the verdict until its name is filled, and goes with its remove button", {tag: ["@demo", "@realizes:u2.records-input", "@realizes:u2.checklist-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Records\" demo page", () => openDemoPage(page, "Records"));
    enter(page, "U2 Demo");
    await session.step(11, "Then Fields records input should have 2 items in total", () => shouldHaveItemsInTotal(page, el("Fields records input"), 2));
    await session.step(12, "And Fields records input should list the items \"Target, Dose\" in this order", () => shouldListInOrder(page, el("Fields records input"), "Target, Dose"));
    await session.step(13, "And value of problem readout should have text \"none\"", () => shouldHaveText(page, el("value of problem readout"), "none"));
    await session.step(14, "When user clicks on \"Add field\" link", () => clickOn(page, el("\"Add field\" link")));
    await session.step(15, "Then Fields records input should have 3 items in total", () => shouldHaveItemsInTotal(page, el("Fields records input"), 3));
    await session.step(16, "And value of problem readout should have text \"Field 3: Name: Value can't be empty\"", () => shouldHaveText(page, el("value of problem readout"), "Field 3: Name: Value can't be empty"));
    await session.step(17, "When user types \"Chemist\" into last text input in Fields records input", () => typeInto(page, "Chemist", el("last text input in Fields records input")));
    await session.step(18, "Then value of problem readout should have text \"none\"", () => shouldHaveText(page, el("value of problem readout"), "none"));
    await session.step(19, "And Fields records input should list the items \"Target, Dose, Chemist\" in this order", () => shouldListInOrder(page, el("Fields records input"), "Target, Dose, Chemist"));
    await session.step(20, "When user clicks on last Remove icon button in Fields records input", () => clickOn(page, el("last Remove icon button in Fields records input")));
    await session.step(21, "Then Fields records input should have 2 items in total", () => shouldHaveItemsInTotal(page, el("Fields records input"), 2));
  });
  test("Stages move with Alt+↑/↓", {tag: ["@demo", "@realizes:u2.records-input", "@realizes:u2.checklist-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Records\" demo page", () => openDemoPage(page, "Records"));
    enter(page, "U2 Demo");
    await session.step(24, "When user presses Alt+ArrowDown in \"Design\" item in Stages records input", () => pressKeyIn(page, "Alt+ArrowDown", el("\"Design\" item in Stages records input")));
    await session.step(25, "Then value of stages readout should have text \"Synthesis, Design, Assay\"", () => shouldHaveText(page, el("value of stages readout"), "Synthesis, Design, Assay"));
    await session.step(26, "And Stages records input should list the items \"Synthesis, Design, Assay\" in this order", () => shouldListInOrder(page, el("Stages records input"), "Synthesis, Design, Assay"));
  });
  test("The checklist's value is its checked items in the order shown", {tag: ["@demo", "@realizes:u2.records-input", "@realizes:u2.checklist-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Records\" demo page", () => openDemoPage(page, "Records"));
    enter(page, "U2 Demo");
    await session.step(29, "Then value of run readout should have text \"Toxicity risks, Descriptors\"", () => shouldHaveText(page, el("value of run readout"), "Toxicity risks, Descriptors"));
    await session.step(30, "And Run checklist should have 4 items in total", () => shouldHaveItemsInTotal(page, el("Run checklist"), 4));
    await session.step(31, "And item \"Toxicity risks\" of Run checklist should be checked", () => itemShouldBe(page, "Toxicity risks", el("Run checklist"), "checked"));
    await session.step(32, "And item \"Solubility\" of Run checklist should be unchecked", () => itemShouldBe(page, "Solubility", el("Run checklist"), "unchecked"));
    await session.step(33, "When user clicks on \"Solubility\" item in Run checklist", () => clickOn(page, el("\"Solubility\" item in Run checklist")));
    await session.step(34, "Then value of run readout should have text \"Toxicity risks, Descriptors, Solubility\"", () => shouldHaveText(page, el("value of run readout"), "Toxicity risks, Descriptors, Solubility"));
    await session.step(35, "When user presses Alt+ArrowUp in \"Solubility\" item in Run checklist", () => pressKeyIn(page, "Alt+ArrowUp", el("\"Solubility\" item in Run checklist")));
    await session.step(36, "Then value of run readout should have text \"Toxicity risks, Solubility, Descriptors\"", () => shouldHaveText(page, el("value of run readout"), "Toxicity risks, Solubility, Descriptors"));
    await session.step(37, "And value of order readout should have text \"Toxicity risks, Solubility, Descriptors, Structural alerts\"", () => shouldHaveText(page, el("value of order readout"), "Toxicity risks, Solubility, Descriptors, Structural alerts"));
  });
});
