/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
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
  const session = feature(test, "features/demo/containers/popups.feature", import.meta.url);
  test("The menu bar writes the status bar", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(11, "When user clicks on File menu item in menu bar", () => clickOn(page, el("File menu item in menu bar")));
    await session.step(12, "And user clicks on \"Save layout\" menu item", () => clickOn(page, el("\"Save layout\" menu item")));
    await session.step(13, "Then status bar should contain text \"Layout saved\"", () => shouldContainText(page, el("status bar"), "Layout saved"));
    await session.step(14, "When user clicks on File menu item in menu bar", () => clickOn(page, el("File menu item in menu bar")));
    await session.step(15, "And user clicks on Export menu item", () => clickOn(page, el("Export menu item")));
    await session.step(16, "And user clicks on \"As CSV\" menu item", () => clickOn(page, el("\"As CSV\" menu item")));
    await session.step(17, "Then status bar should contain text \"Exported as CSV\"", () => shouldContainText(page, el("status bar"), "Exported as CSV"));
    await session.step(18, "When user clicks on View menu item in menu bar", () => clickOn(page, el("View menu item in menu bar")));
    await session.step(19, "And user clicks on Auto-refresh menu item", () => clickOn(page, el("Auto-refresh menu item")));
    await session.step(20, "Then status bar should contain text \"Auto-refresh on\"", () => shouldContainText(page, el("status bar"), "Auto-refresh on"));
  });
  test("A popup menu with a shortcut hint and a disabled item", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(23, "When user clicks on \"Open menu\" button", () => clickOn(page, el("\"Open menu\" button")));
    await session.step(24, "Then menu should be visible", () => shouldBe(page, el("menu"), "visible"));
    await session.step(25, "And \"Add to favorites\" menu item should contain text \"Ctrl+D\"", () => shouldContainText(page, el("\"Add to favorites\" menu item"), "Ctrl+D"));
    await session.step(26, "And Delete menu item should be disabled", () => shouldBe(page, el("Delete menu item"), "disabled"));
    await session.step(27, "When user clicks on Rename menu item", () => clickOn(page, el("Rename menu item")));
    await session.step(28, "Then status bar should contain text \"Renamed\"", () => shouldContainText(page, el("status bar"), "Renamed"));
    await session.step(29, "And menu should be hidden", () => shouldBe(page, el("menu"), "hidden"));
  });
  test("A context menu on a panel", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(32, "When user right-clicks on \"Right-click anywhere in this panel for a context menu.\" text", () => rightClickOn(page, el("\"Right-click anywhere in this panel for a context menu.\" text")));
    await session.step(33, "And user clicks on \"Select all\" menu item", () => clickOn(page, el("\"Select all\" menu item")));
    await session.step(34, "Then status bar should contain text \"Selected all\"", () => shouldContainText(page, el("status bar"), "Selected all"));
  });
  test("A tooltip evaluated at show time", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(37, "When user hovers over \"Hover me for a tooltip.\" text", () => hoverOver(page, el("\"Hover me for a tooltip.\" text")));
    await session.step(38, "Then tooltip should contain text \"One shared tooltip element, evaluated at show time\"", () => shouldContainText(page, el("tooltip"), "One shared tooltip element, evaluated at show time"));
  });
  test("The dialog accepts and cancels", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(41, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(42, "Then dialog should be visible", () => shouldBe(page, el("dialog"), "visible"));
    await session.step(43, "And title of dialog should have text \"u2 dialog\"", () => shouldHaveText(page, el("title of dialog"), "u2 dialog"));
    await session.step(44, "And Name input in dialog should have value \"u2\"", () => shouldHaveValue(page, el("Name input in dialog"), "u2"));
    await session.step(45, "When user types \"bdd\" into Name input in dialog", () => typeInto(page, "bdd", el("Name input in dialog")));
    await session.step(46, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await session.step(47, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
    await session.step(48, "And status bar should contain text \"Dialog OK, name = bdd\"", () => shouldContainText(page, el("status bar"), "Dialog OK, name = bdd"));
    await session.step(49, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(50, "And user closes dialog", () => close(page, el("dialog")));
    await session.step(51, "Then status bar should contain text \"Dialog cancelled\"", () => shouldContainText(page, el("status bar"), "Dialog cancelled"));
  });
});
