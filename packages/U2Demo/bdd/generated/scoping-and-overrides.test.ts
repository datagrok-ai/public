/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scoping-and-overrides.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/demo.js';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openWorkbench} from '../bindings/steps.js';
import {clickOn, close, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Scoping and overrides", () => {
  const session = feature(test, "features/scoping-and-overrides.feature", import.meta.url);
  test("The same button name in two scopes", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(12, "When user clicks on save button in alignment panel", () => clickOn(page, el("save button in alignment panel")));
    await session.step(13, "Then notification should contain text \"Saved from the form\"", () => shouldContainText(page, el("notification"), "Saved from the form"));
    await session.step(14, "When user clicks on save button in toolbar", () => clickOn(page, el("save button in toolbar")));
    await session.step(15, "Then notification should contain text \"Saved from the toolbar\"", () => shouldContainText(page, el("notification"), "Saved from the toolbar"));
  });
  test("A registered whole phrase overrides composition", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(18, "When user clicks on save button inside toolbar", () => clickOn(page, el("save button inside toolbar")));
    await session.step(19, "Then notification should contain text \"Saved from the toolbar\"", () => shouldContainText(page, el("notification"), "Saved from the toolbar"));
  });
  test("Generic kinds inside a registered scope", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(22, "When user clicks on run button in toolbar", () => clickOn(page, el("run button in toolbar")));
    await session.step(23, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(24, "When user closes MSA dialog", () => close(page, el("MSA dialog")));
    await session.step(25, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
  });
  test("Ordinals among the matches", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(28, "When user clicks on run button in toolbar", () => clickOn(page, el("run button in toolbar")));
    await session.step(29, "And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await session.step(30, "Then second item in aligned sequences list should contain text \"2\"", () => shouldContainText(page, el("second item in aligned sequences list"), "2"));
    await session.step(31, "And last item in aligned sequences list should contain text \"5\"", () => shouldContainText(page, el("last item in aligned sequences list"), "5"));
  });
  test("A toggle reveals a panel", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(34, "Then settings panel should be hidden", () => shouldBe(page, el("settings panel"), "hidden"));
    await session.step(35, "When user clicks on settings button in toolbar", () => clickOn(page, el("settings button in toolbar")));
    await session.step(36, "Then settings panel should be visible", () => shouldBe(page, el("settings panel"), "visible"));
    await session.step(37, "And theme input in settings panel should have value \"light\"", () => shouldHaveValue(page, el("theme input in settings panel"), "light"));
  });
  test("Platform names keep their platform meaning inside the workbench", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(40, "Then toolbar should be visible", () => shouldBe(page, el("toolbar"), "visible"));
    await session.step(41, "And MSA workbench should be visible", () => shouldBe(page, el("MSA workbench"), "visible"));
    await session.step(42, "And toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(43, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
  });
});
