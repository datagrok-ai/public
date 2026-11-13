/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/table-chrome.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.ribbon-group]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {windowSize} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {clickOn, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal, shouldReportReading} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("TableView chrome", () => {
  const session = feature(test, "features/demo/platform/table-chrome.feature", import.meta.url);
  test("The group, its icons and buttons; a dialog pinned to the view", {tag: ["@demo", "@realizes:u2.ribbon-group"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"TableView chrome\" demo page", () => openDemoPage(page, "TableView chrome"));
    enter(page, "U2 Demo");
    await session.step(9, "When user clicks on \"Open the campaign view\" button", () => clickOn(page, el("\"Open the campaign view\" button")));
    await session.step(10, "Then the \"KHL-3\" view should be current", () => viewIsCurrent(page, "KHL-3"));
    await session.step(13, "Then \"demo-campaign-actions\" ribbon group should have 9 items in total", () => shouldHaveItemsInTotal(page, el("\"demo-campaign-actions\" ribbon group"), 9));
    await session.step(14, "And \"demo-campaign-actions\" ribbon group should report \"folded\" as \"false\"", () => shouldReportReading(page, el("\"demo-campaign-actions\" ribbon group"), "folded", "false"));
    await session.step(15, "When user clicks on \"Add new row\" icon button", () => clickOn(page, el("\"Add new row\" icon button")));
    await session.step(16, "Then status bar should contain text \"of 4\"", () => shouldContainText(page, el("status bar"), "of 4"));
    await session.step(17, "When user clicks on \"Progress tracker\" button", () => clickOn(page, el("\"Progress tracker\" button")));
    await session.step(18, "Then \"Progress tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress tracker\" dialog"), "visible"));
    await session.step(19, "When user switches to the \"U2 Demo\" view", () => switchView(page, "U2 Demo"));
    await session.step(20, "Then \"Progress tracker\" dialog should be hidden", () => shouldBe(page, el("\"Progress tracker\" dialog"), "hidden"));
    await session.step(21, "When user switches to the \"KHL-3\" view", () => switchView(page, "KHL-3"));
    await session.step(22, "Then \"Progress tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress tracker\" dialog"), "visible"));
    await session.step(23, "When user clicks on CANCEL button in \"Progress tracker\" dialog", () => clickOn(page, el("CANCEL button in \"Progress tracker\" dialog")));
    await session.step(24, "Then \"Progress tracker\" dialog should be hidden", () => shouldBe(page, el("\"Progress tracker\" dialog"), "hidden"));
  });
  test("Narrow, the icons are one \"…\" menu and the buttons stay", {tag: ["@demo", "@realizes:u2.ribbon-group"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"TableView chrome\" demo page", () => openDemoPage(page, "TableView chrome"));
    enter(page, "U2 Demo");
    await session.step(9, "When user clicks on \"Open the campaign view\" button", () => clickOn(page, el("\"Open the campaign view\" button")));
    await session.step(10, "Then the \"KHL-3\" view should be current", () => viewIsCurrent(page, "KHL-3"));
    await session.step(27, "Given the window is 1024 by 768 px", () => windowSize(page, 1024, 768));
    await session.step(28, "Then \"demo-campaign-actions\" ribbon group should report \"folded\" as \"true\"", () => shouldReportReading(page, el("\"demo-campaign-actions\" ribbon group"), "folded", "true"));
    await session.step(29, "And \"Submit...\" button should be visible", () => shouldBe(page, el("\"Submit...\" button"), "visible"));
    await session.step(30, "And \"Progress tracker\" button should be visible", () => shouldBe(page, el("\"Progress tracker\" button"), "visible"));
    await session.step(31, "When user clicks on \"More actions\" icon button", () => clickOn(page, el("\"More actions\" icon button")));
    await session.step(32, "And user clicks on \"Merge another table into the campaign...\" menu item", () => clickOn(page, el("\"Merge another table into the campaign...\" menu item")));
    await session.step(33, "Then status bar should contain text \"Merge\"", () => shouldContainText(page, el("status bar"), "Merge"));
  });
  test("The crumb goes back and closes the campaign view", {tag: ["@demo", "@realizes:u2.ribbon-group"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"TableView chrome\" demo page", () => openDemoPage(page, "TableView chrome"));
    enter(page, "U2 Demo");
    await session.step(9, "When user clicks on \"Open the campaign view\" button", () => clickOn(page, el("\"Open the campaign view\" button")));
    await session.step(10, "Then the \"KHL-3\" view should be current", () => viewIsCurrent(page, "KHL-3"));
    await session.step(36, "When user clicks on \"U2 Demo\" breadcrumb", () => clickOn(page, el("\"U2 Demo\" breadcrumb")));
    await session.step(37, "Then the \"U2 Demo\" view should be current", () => viewIsCurrent(page, "U2 Demo"));
    await session.step(38, "And value of campaign view readout should have text \"closed\"", () => shouldHaveText(page, el("value of campaign view readout"), "closed"));
  });
});
