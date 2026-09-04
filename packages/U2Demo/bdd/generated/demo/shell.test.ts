/* ---
generated: features/demo/shell.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.app-view]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../bindings/demo.js';
import {clickOn, collapse, expand, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The demo shell", () => {
  const session = feature(test);
  test("Navigating through the tree and the links", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await test.step("Then status bar should contain text \"Start / Overview\"", () => shouldContainText(page, el("status bar"), "Start / Overview"));
    await test.step("When user clicks on \"Forms / Form\" link", () => clickOn(page, el("\"Forms / Form\" link")));
    await test.step("Then status bar should contain text \"Forms / Form\"", () => shouldContainText(page, el("status bar"), "Forms / Form"));
    await test.step("And Form tree node in demo navigation should be selected", () => shouldBe(page, el("Form tree node in demo navigation"), "selected"));
    await test.step("When user clicks on Cards tree node in demo navigation", () => clickOn(page, el("Cards tree node in demo navigation")));
    await test.step("Then status bar should contain text \"Display / Cards\"", () => shouldContainText(page, el("status bar"), "Display / Cards"));
    await test.step("And Revenue stat card should be visible", () => shouldBe(page, el("Revenue stat card"), "visible"));
    await test.step("When user collapses Display tree node in demo navigation", () => collapse(page, el("Display tree node in demo navigation")));
    await test.step("Then Cards tree node in demo navigation should be hidden", () => shouldBe(page, el("Cards tree node in demo navigation"), "hidden"));
    await test.step("When user expands Display tree node in demo navigation", () => expand(page, el("Display tree node in demo navigation")));
    await test.step("Then Cards tree node in demo navigation should be visible", () => shouldBe(page, el("Cards tree node in demo navigation"), "visible"));
  });
  test("The ribbon's demo tools open and close tables", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"Demo tools\" dropdown button", () => clickOn(page, el("\"Demo tools\" dropdown button")));
    await test.step("And user clicks on \"Add demog table\" menu item", () => clickOn(page, el("\"Add demog table\" menu item")));
    await test.step("Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await test.step("When user clicks on \"U2 Demo\" view", () => clickOn(page, el("\"U2 Demo\" view")));
    await test.step("And user clicks on \"Demo tools\" dropdown button", () => clickOn(page, el("\"Demo tools\" dropdown button")));
    await test.step("And user clicks on \"Close demo tables\" menu item", () => clickOn(page, el("\"Close demo tables\" menu item")));
    await test.step("Then grid should be hidden", () => shouldBe(page, el("grid"), "hidden"));
  });
});
