/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
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
  const session = feature(test, "features/demo/shell.feature", import.meta.url);
  test("Navigating through the tree and the links", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await session.step(8, "Then status bar should contain text \"Start / Overview\"", () => shouldContainText(page, el("status bar"), "Start / Overview"));
    await session.step(9, "When user clicks on \"Forms / Form\" link", () => clickOn(page, el("\"Forms / Form\" link")));
    await session.step(10, "Then status bar should contain text \"Forms / Form\"", () => shouldContainText(page, el("status bar"), "Forms / Form"));
    await session.step(11, "And Form tree node in demo navigation should be selected", () => shouldBe(page, el("Form tree node in demo navigation"), "selected"));
    await session.step(12, "When user clicks on Cards tree node in demo navigation", () => clickOn(page, el("Cards tree node in demo navigation")));
    await session.step(13, "Then status bar should contain text \"Display / Cards\"", () => shouldContainText(page, el("status bar"), "Display / Cards"));
    await session.step(14, "And Revenue stat card should be visible", () => shouldBe(page, el("Revenue stat card"), "visible"));
    await session.step(15, "When user collapses Display tree node in demo navigation", () => collapse(page, el("Display tree node in demo navigation")));
    await session.step(16, "Then Cards tree node in demo navigation should be hidden", () => shouldBe(page, el("Cards tree node in demo navigation"), "hidden"));
    await session.step(17, "When user expands Display tree node in demo navigation", () => expand(page, el("Display tree node in demo navigation")));
    await session.step(18, "Then Cards tree node in demo navigation should be visible", () => shouldBe(page, el("Cards tree node in demo navigation"), "visible"));
  });
  test("The ribbon's demo tools open and close tables", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await session.step(22, "When user clicks on \"Demo tools\" dropdown button", () => clickOn(page, el("\"Demo tools\" dropdown button")));
    await session.step(23, "And user clicks on \"Add demog table\" menu item", () => clickOn(page, el("\"Add demog table\" menu item")));
    await session.step(24, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(25, "When user clicks on \"U2 Demo\" view", () => clickOn(page, el("\"U2 Demo\" view")));
    await session.step(26, "And user clicks on \"Demo tools\" dropdown button", () => clickOn(page, el("\"Demo tools\" dropdown button")));
    await session.step(27, "And user clicks on \"Close demo tables\" menu item", () => clickOn(page, el("\"Close demo tables\" menu item")));
    await session.step(28, "Then grid should be hidden", () => shouldBe(page, el("grid"), "hidden"));
  });
});
