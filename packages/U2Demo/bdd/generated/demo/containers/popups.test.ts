/* ---
generated: features/demo/containers/popups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.menu, u2.menu-bar, u2.dialog, u2.tooltip]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, close, hoverOver, rightClickOn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Popups", () => {
  const session = feature(test);
  test("The menu bar writes the status bar", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on File menu item in menu bar", () => clickOn(page, el("File menu item in menu bar")));
    await test.step("And user clicks on \"Save layout\" menu item", () => clickOn(page, el("\"Save layout\" menu item")));
    await test.step("Then status bar should contain text \"Layout saved\"", () => shouldContainText(page, el("status bar"), "Layout saved"));
    await test.step("When user clicks on File menu item in menu bar", () => clickOn(page, el("File menu item in menu bar")));
    await test.step("And user clicks on Export menu item", () => clickOn(page, el("Export menu item")));
    await test.step("And user clicks on \"As CSV\" menu item", () => clickOn(page, el("\"As CSV\" menu item")));
    await test.step("Then status bar should contain text \"Exported as CSV\"", () => shouldContainText(page, el("status bar"), "Exported as CSV"));
    await test.step("When user clicks on View menu item in menu bar", () => clickOn(page, el("View menu item in menu bar")));
    await test.step("And user clicks on Auto-refresh menu item", () => clickOn(page, el("Auto-refresh menu item")));
    await test.step("Then status bar should contain text \"Auto-refresh on\"", () => shouldContainText(page, el("status bar"), "Auto-refresh on"));
  });
  test("A popup menu with a shortcut hint and a disabled item", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"Open menu\" button", () => clickOn(page, el("\"Open menu\" button")));
    await test.step("Then menu should be visible", () => shouldBe(page, el("menu"), "visible"));
    await test.step("And \"Add to favorites\" menu item should contain text \"Ctrl+D\"", () => shouldContainText(page, el("\"Add to favorites\" menu item"), "Ctrl+D"));
    await test.step("And Delete menu item should be disabled", () => shouldBe(page, el("Delete menu item"), "disabled"));
    await test.step("When user clicks on Rename menu item", () => clickOn(page, el("Rename menu item")));
    await test.step("Then status bar should contain text \"Renamed\"", () => shouldContainText(page, el("status bar"), "Renamed"));
    await test.step("And menu should be hidden", () => shouldBe(page, el("menu"), "hidden"));
  });
  test("A context menu on a panel", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await test.step("When user right-clicks on \"Right-click anywhere in this panel for a context menu.\" text", () => rightClickOn(page, el("\"Right-click anywhere in this panel for a context menu.\" text")));
    await test.step("And user clicks on \"Select all\" menu item", () => clickOn(page, el("\"Select all\" menu item")));
    await test.step("Then status bar should contain text \"Selected all\"", () => shouldContainText(page, el("status bar"), "Selected all"));
  });
  test("A tooltip evaluated at show time", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await test.step("When user hovers over \"Hover me for a tooltip.\" text", () => hoverOver(page, el("\"Hover me for a tooltip.\" text")));
    await test.step("Then tooltip should contain text \"One shared tooltip element, evaluated at show time\"", () => shouldContainText(page, el("tooltip"), "One shared tooltip element, evaluated at show time"));
  });
  test("The dialog accepts and cancels", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await test.step("When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await test.step("Then dialog should be visible", () => shouldBe(page, el("dialog"), "visible"));
    await test.step("And title of dialog should have text \"u2 dialog\"", () => shouldHaveText(page, el("title of dialog"), "u2 dialog"));
    await test.step("And Name input in dialog should have value \"u2\"", () => shouldHaveValue(page, el("Name input in dialog"), "u2"));
    await test.step("When user types \"bdd\" into Name input in dialog", () => typeInto(page, "bdd", el("Name input in dialog")));
    await test.step("And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await test.step("Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
    await test.step("And status bar should contain text \"Dialog OK, name = bdd\"", () => shouldContainText(page, el("status bar"), "Dialog OK, name = bdd"));
    await test.step("When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await test.step("And user closes dialog", () => close(page, el("dialog")));
    await test.step("Then status bar should contain text \"Dialog cancelled\"", () => shouldContainText(page, el("status bar"), "Dialog cancelled"));
  });
});
