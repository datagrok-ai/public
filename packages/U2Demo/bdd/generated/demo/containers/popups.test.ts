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
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {atLeastWidePx, insideWindow, tallPx, widePx, windowSize} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {clickOn, close, dragBy, hoverOver, pressKey, rightClickOn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
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
    await session.step(27, "And \"Add to favorites\" menu item should be 24 px tall", () => tallPx(page, el("\"Add to favorites\" menu item"), 24));
    await session.step(28, "And menu should be at least 160 px wide", () => atLeastWidePx(page, el("menu"), 160));
    await session.step(29, "When user clicks on Rename menu item", () => clickOn(page, el("Rename menu item")));
    await session.step(30, "Then status bar should contain text \"Renamed\"", () => shouldContainText(page, el("status bar"), "Renamed"));
    await session.step(31, "And menu should be hidden", () => shouldBe(page, el("menu"), "hidden"));
  });
  test("A context menu on a panel", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(34, "When user right-clicks on \"Right-click anywhere in this panel for a context menu.\" text", () => rightClickOn(page, el("\"Right-click anywhere in this panel for a context menu.\" text")));
    await session.step(35, "And user clicks on \"Select all\" menu item", () => clickOn(page, el("\"Select all\" menu item")));
    await session.step(36, "Then status bar should contain text \"Selected all\"", () => shouldContainText(page, el("status bar"), "Selected all"));
  });
  test("A tooltip evaluated at show time", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(39, "When user hovers over \"Hover me for a tooltip.\" text", () => hoverOver(page, el("\"Hover me for a tooltip.\" text")));
    await session.step(40, "Then tooltip should contain text \"One shared tooltip element, evaluated at show time\"", () => shouldContainText(page, el("tooltip"), "One shared tooltip element, evaluated at show time"));
  });
  test("The dialog accepts and cancels", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(43, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(44, "Then dialog should be visible", () => shouldBe(page, el("dialog"), "visible"));
    await session.step(45, "And title of dialog should have text \"u2 dialog\"", () => shouldHaveText(page, el("title of dialog"), "u2 dialog"));
    await session.step(46, "And Name input in dialog should have value \"u2\"", () => shouldHaveValue(page, el("Name input in dialog"), "u2"));
    await session.step(47, "And editor of Opacity input in dialog should be at least 200 px wide", () => atLeastWidePx(page, el("editor of Opacity input in dialog"), 200));
    await session.step(48, "And editor of Targets input in dialog should be at least 150 px wide", () => atLeastWidePx(page, el("editor of Targets input in dialog"), 150));
    await session.step(49, "When user types \"bdd\" into Name input in dialog", () => typeInto(page, "bdd", el("Name input in dialog")));
    await session.step(50, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await session.step(51, "Then dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
    await session.step(52, "And status bar should contain text \"Dialog OK, name = bdd\"", () => shouldContainText(page, el("status bar"), "Dialog OK, name = bdd"));
    await session.step(53, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(54, "And user closes dialog", () => close(page, el("dialog")));
    await session.step(55, "Then status bar should contain text \"Dialog cancelled\"", () => shouldContainText(page, el("status bar"), "Dialog cancelled"));
  });
  test("The dialog chrome reads as the platform's", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(58, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(59, "Then footer of dialog should be 45 px tall", () => tallPx(page, el("footer of dialog"), 45));
    await session.step(60, "And close button of dialog should be 18 px wide", () => widePx(page, el("close button of dialog"), 18));
    await session.step(61, "When user clicks on Help button in dialog", () => clickOn(page, el("Help button in dialog")));
    await session.step(62, "Then status bar should contain text \"Dialog help\"", () => shouldContainText(page, el("status bar"), "Dialog help"));
    await session.step(63, "And dialog should be visible", () => shouldBe(page, el("dialog"), "visible"));
  });
  test("The dialog selects the first editor's text on show", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(66, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(67, "And user presses 7", () => pressKey(page, "7"));
    await session.step(68, "Then Name input in dialog should have value \"7\"", () => shouldHaveValue(page, el("Name input in dialog"), "7"));
  });
  test("A dialog whose content grows stays on screen", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(71, "Given the window is 1024 by 500 px", () => windowSize(page, 1024, 500));
    await session.step(72, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(73, "And user clicks on \"Add rows\" button in dialog", () => clickOn(page, el("\"Add rows\" button in dialog")));
    await session.step(74, "Then footer of dialog should be inside the window", () => insideWindow(page, el("footer of dialog")));
    await session.step(75, "And OK button in dialog should be inside the window", () => insideWindow(page, el("OK button in dialog")));
  });
  test("A dialog dragged right keeps its size when the window narrows", {tag: ["@demo", "@realizes:u2.menu", "@realizes:u2.menu-bar", "@realizes:u2.dialog", "@realizes:u2.tooltip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Popups\" demo page", () => openDemoPage(page, "Popups"));
    enter(page, "U2 Demo");
    await session.step(78, "When user clicks on \"Open dialog\" button", () => clickOn(page, el("\"Open dialog\" button")));
    await session.step(79, "And user drags title of dialog by 600 and 0 pixels", () => dragBy(page, el("title of dialog"), 600, 0));
    await session.step(80, "And the window is 900 by 900 px", () => windowSize(page, 900, 900));
    await session.step(81, "Then dialog should be at least 360 px wide", () => atLeastWidePx(page, el("dialog"), 360));
    await session.step(82, "And dialog should be inside the window", () => insideWindow(page, el("dialog")));
  });
});
