/* ---
generated: features/demo/display/cards.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.card, u2.stat-card]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, shouldBe, shouldContainText, shouldHaveText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Cards and stat cards", () => {
  const session = feature(test);
  test("A clickable card counts its clicks", {tag: ["@demo", "@realizes:u2.card", "@realizes:u2.stat-card"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Cards\" demo page", () => openDemoPage(page, "Cards"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"Open the compound\" card", () => clickOn(page, el("\"Open the compound\" card")));
    await test.step("And user clicks on \"Open the compound\" card", () => clickOn(page, el("\"Open the compound\" card")));
    await test.step("Then value of clicks readout should have text \"2\"", () => shouldHaveText(page, el("value of clicks readout"), "2"));
  });
  test("Selectable cards adopt the signals they are given", {tag: ["@demo", "@realizes:u2.card", "@realizes:u2.stat-card"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Cards\" demo page", () => openDemoPage(page, "Cards"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on Caffeine card", () => clickOn(page, el("Caffeine card")));
    await test.step("Then Caffeine card should be selected", () => shouldBe(page, el("Caffeine card"), "selected"));
    await test.step("And value of selected readout should have text \"Caffeine\"", () => shouldHaveText(page, el("value of selected readout"), "Caffeine"));
    await test.step("When user clicks on second Aspirin card", () => clickOn(page, el("second Aspirin card")));
    await test.step("Then value of selected readout should have text \"Aspirin, Caffeine\"", () => shouldHaveText(page, el("value of selected readout"), "Aspirin, Caffeine"));
    await test.step("When user clicks on \"Select all\" button", () => clickOn(page, el("\"Select all\" button")));
    await test.step("Then value of selected readout should have text \"Aspirin, Caffeine, Ibuprofen\"", () => shouldHaveText(page, el("value of selected readout"), "Aspirin, Caffeine, Ibuprofen"));
    await test.step("When user clicks on Clear button", () => clickOn(page, el("Clear button")));
    await test.step("Then value of selected readout should have text \"(none)\"", () => shouldHaveText(page, el("value of selected readout"), "(none)"));
    await test.step("And Caffeine card should not be selected", () => shouldNotBe(page, el("Caffeine card"), "selected"));
  });
  test("Stat cards follow their signals", {tag: ["@demo", "@realizes:u2.card", "@realizes:u2.stat-card"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Cards\" demo page", () => openDemoPage(page, "Cards"));
    enter(page, "U2 Demo");
    await test.step("Then Revenue stat card should contain text \"1.24M\"", () => shouldContainText(page, el("Revenue stat card"), "1.24M"));
    await test.step("When user clicks on \"+120k\" button", () => clickOn(page, el("\"+120k\" button")));
    await test.step("Then Revenue stat card should contain text \"1.36M\"", () => shouldContainText(page, el("Revenue stat card"), "1.36M"));
    await test.step("And value of delta readout should have text \"0.12\"", () => shouldHaveText(page, el("value of delta readout"), "0.12"));
    await test.step("When user clicks on \"-120k\" button", () => clickOn(page, el("\"-120k\" button")));
    await test.step("Then Revenue stat card should contain text \"1.24M\"", () => shouldContainText(page, el("Revenue stat card"), "1.24M"));
    await test.step("And value of delta readout should have text \"-0.09\"", () => shouldHaveText(page, el("value of delta readout"), "-0.09"));
  });
});
