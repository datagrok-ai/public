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
  const session = feature(test);
  test("The same button name in two scopes", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user clicks on save button in alignment panel", () => clickOn(page, el("save button in alignment panel")));
    await test.step("Then notification should contain text \"Saved from the form\"", () => shouldContainText(page, el("notification"), "Saved from the form"));
    await test.step("When user clicks on save button in toolbar", () => clickOn(page, el("save button in toolbar")));
    await test.step("Then notification should contain text \"Saved from the toolbar\"", () => shouldContainText(page, el("notification"), "Saved from the toolbar"));
  });
  test("A registered whole phrase overrides composition", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user clicks on save button inside toolbar", () => clickOn(page, el("save button inside toolbar")));
    await test.step("Then notification should contain text \"Saved from the toolbar\"", () => shouldContainText(page, el("notification"), "Saved from the toolbar"));
  });
  test("Generic kinds inside a registered scope", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user clicks on run button in toolbar", () => clickOn(page, el("run button in toolbar")));
    await test.step("Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await test.step("When user closes MSA dialog", () => close(page, el("MSA dialog")));
    await test.step("Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
  });
  test("Ordinals among the matches", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user clicks on run button in toolbar", () => clickOn(page, el("run button in toolbar")));
    await test.step("And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await test.step("Then second item in aligned sequences list should contain text \"2\"", () => shouldContainText(page, el("second item in aligned sequences list"), "2"));
    await test.step("And last item in aligned sequences list should contain text \"5\"", () => shouldContainText(page, el("last item in aligned sequences list"), "5"));
  });
  test("A toggle reveals a panel", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("Then settings panel should be hidden", () => shouldBe(page, el("settings panel"), "hidden"));
    await test.step("When user clicks on settings button in toolbar", () => clickOn(page, el("settings button in toolbar")));
    await test.step("Then settings panel should be visible", () => shouldBe(page, el("settings panel"), "visible"));
    await test.step("And theme input in settings panel should have value \"light\"", () => shouldHaveValue(page, el("theme input in settings panel"), "light"));
  });
  test("Platform names keep their platform meaning inside the workbench", {tag: ["@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("Then toolbar should be visible", () => shouldBe(page, el("toolbar"), "visible"));
    await test.step("And MSA workbench should be visible", () => shouldBe(page, el("MSA workbench"), "visible"));
    await test.step("And toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await test.step("And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
  });
});
