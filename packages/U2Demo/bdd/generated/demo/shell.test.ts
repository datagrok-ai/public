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
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {navigationShown, openDemoPage} from '../../bindings/demo.js';
import {atLeastWidePx, leftOf, widePx, windowSize} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {clickOn, collapse, expand, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, switchView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The demo shell", () => {
  const session = feature(test, "features/demo/shell.feature", import.meta.url);
  test("Navigating through the tree and the links", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await session.step(8, "And the demo navigation is shown", () => navigationShown(page));
    await session.step(9, "Then demo navigation should be at least 148 px wide", () => atLeastWidePx(page, el("demo navigation"), 148));
    await session.step(10, "And demo navigation should be left of demo page", () => leftOf(page, el("demo navigation"), el("demo page")));
    await session.step(11, "And demo page should be at least 280 px wide", () => atLeastWidePx(page, el("demo page"), 280));
    await session.step(12, "And status bar should contain text \"Start / Overview\"", () => shouldContainText(page, el("status bar"), "Start / Overview"));
    await session.step(13, "When user clicks on \"Forms / Form\" link", () => clickOn(page, el("\"Forms / Form\" link")));
    await session.step(14, "Then status bar should contain text \"Forms / Form\"", () => shouldContainText(page, el("status bar"), "Forms / Form"));
    await session.step(15, "And Form tree node in demo navigation should be selected", () => shouldBe(page, el("Form tree node in demo navigation"), "selected"));
    await session.step(16, "When user clicks on Cards tree node in demo navigation", () => clickOn(page, el("Cards tree node in demo navigation")));
    await session.step(17, "Then status bar should contain text \"Display / Cards\"", () => shouldContainText(page, el("status bar"), "Display / Cards"));
    await session.step(18, "And Revenue stat card should be visible", () => shouldBe(page, el("Revenue stat card"), "visible"));
    await session.step(19, "When user collapses Display tree node in demo navigation", () => collapse(page, el("Display tree node in demo navigation")));
    await session.step(20, "Then Cards tree node in demo navigation should be hidden", () => shouldBe(page, el("Cards tree node in demo navigation"), "hidden"));
    await session.step(21, "When user expands Display tree node in demo navigation", () => expand(page, el("Display tree node in demo navigation")));
    await session.step(22, "Then Cards tree node in demo navigation should be visible", () => shouldBe(page, el("Cards tree node in demo navigation"), "visible"));
  });
  test("A pane too narrow for the nav folds it away, and the ribbon brings the rail back", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the window is 1024 by 768 px", () => windowSize(page, 1024, 768));
    await session.step(26, "And user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await session.step(27, "Then demo navigation should be hidden", () => shouldBe(page, el("demo navigation"), "hidden"));
    await session.step(28, "When user clicks on \"Show navigation\" icon button", () => clickOn(page, el("\"Show navigation\" icon button")));
    await session.step(29, "Then demo navigation should be visible", () => shouldBe(page, el("demo navigation"), "visible"));
    await session.step(30, "And demo navigation should be 148 px wide", () => widePx(page, el("demo navigation"), 148));
    await session.step(31, "And demo navigation should be left of demo page", () => leftOf(page, el("demo navigation"), el("demo page")));
  });
  test("On a wide screen the nav keeps its width in pixels", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given the window is 1920 by 1080 px", () => windowSize(page, 1920, 1080));
    await session.step(35, "And user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await session.step(36, "Then demo navigation should be 180 px wide", () => widePx(page, el("demo navigation"), 180));
  });
  test("The context panel shows the source of the page, whatever file of the demo it lives in [page=Records, factory=recordsPage]", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(39, "Given user opens the \"Records\" demo page", () => openDemoPage(page, "Records"));
    enter(page, "U2 Demo");
    await session.step(40, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(41, "Then context panel should contain text \"export function recordsPage(\"", () => shouldContainText(page, el("context panel"), "export function recordsPage("));
    await session.step(42, "And context panel should not contain text \"Not bundled\"", () => shouldNotContainText(page, el("context panel"), "Not bundled"));
  });
  test("The context panel shows the source of the page, whatever file of the demo it lives in [page=TableView chrome, factory=tableChromePage]", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(39, "Given user opens the \"TableView chrome\" demo page", () => openDemoPage(page, "TableView chrome"));
    enter(page, "U2 Demo");
    await session.step(40, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(41, "Then context panel should contain text \"export function tableChromePage(\"", () => shouldContainText(page, el("context panel"), "export function tableChromePage("));
    await session.step(42, "And context panel should not contain text \"Not bundled\"", () => shouldNotContainText(page, el("context panel"), "Not bundled"));
  });
  test("The context panel shows the source of the page, whatever file of the demo it lives in [page=Tables, factory=tablesPage]", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(39, "Given user opens the \"Tables\" demo page", () => openDemoPage(page, "Tables"));
    enter(page, "U2 Demo");
    await session.step(40, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(41, "Then context panel should contain text \"export function tablesPage(\"", () => shouldContainText(page, el("context panel"), "export function tablesPage("));
    await session.step(42, "And context panel should not contain text \"Not bundled\"", () => shouldNotContainText(page, el("context panel"), "Not bundled"));
  });
  test("The ribbon's demo tools open and close tables", {tag: ["@demo", "@realizes:u2.app-view"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(51, "Given user opens the \"Overview\" demo page", () => openDemoPage(page, "Overview"));
    enter(page, "U2 Demo");
    await session.step(52, "When user clicks on \"Demo tools\" dropdown button", () => clickOn(page, el("\"Demo tools\" dropdown button")));
    await session.step(53, "And user clicks on \"Add demog table\" menu item", () => clickOn(page, el("\"Add demog table\" menu item")));
    await session.step(54, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(55, "When user switches to the \"U2 Demo\" view", () => switchView(page, "U2 Demo"));
    await session.step(56, "And user clicks on \"Demo tools\" dropdown button", () => clickOn(page, el("\"Demo tools\" dropdown button")));
    await session.step(57, "And user clicks on \"Close demo tables\" menu item", () => clickOn(page, el("\"Close demo tables\" menu item")));
    await session.step(58, "Then grid should be hidden", () => shouldBe(page, el("grid"), "hidden"));
  });
});
