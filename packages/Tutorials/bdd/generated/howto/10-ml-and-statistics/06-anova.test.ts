/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/06-anova.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/06-anova.feature", import.meta.url);
  test("Compare groups with ANOVA", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When the video says \"ANOVA checks whether age differs between six disease groups\" at the bottom", () => videoSaysAt(page, "ANOVA checks whether age differs between six disease groups", "bottom"));
    await session.step(22, "And user picks \"ML > Analyze > Group Comparison > ANOVA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > ANOVA..."), undefined, "ANOVA is under ML > Analyze > Group Comparison");
    await session.step(23, "Then \"ANOVA\" dialog should be visible", () => shouldBe(page, el("\"ANOVA\" dialog"), "visible"));
    await session.step(25, "When user selects \"DIS_POP\" in Category input in \"ANOVA\" dialog", () => selectIn(page, "DIS_POP", el("Category input in \"ANOVA\" dialog")), undefined, "Category: the column that splits rows into groups");
    await session.step(26, "Then editor of Category input in \"ANOVA\" dialog should have text \"DIS_POP\"", () => shouldHaveText(page, el("editor of Category input in \"ANOVA\" dialog"), "DIS_POP"));
    await session.step(28, "When user hovers over editor of Feature input in \"ANOVA\" dialog", () => hoverOver(page, el("editor of Feature input in \"ANOVA\" dialog")), undefined, "Feature: the number to compare between the groups");
    await session.step(29, "Then editor of Feature input in \"ANOVA\" dialog should have text \"AGE\"", () => shouldHaveText(page, el("editor of Feature input in \"ANOVA\" dialog"), "AGE"));
    await session.step(31, "When user hovers over Alpha input in \"ANOVA\" dialog", () => hoverOver(page, el("Alpha input in \"ANOVA\" dialog")), undefined, "Alpha: how sure you want to be, 0.05 is usual");
    await session.step(32, "Then Alpha input in \"ANOVA\" dialog should have value \"0.05\"", () => shouldHaveValue(page, el("Alpha input in \"ANOVA\" dialog"), "0.05"));
    await session.step(33, "When user clicks on Run button in \"ANOVA\" dialog", () => clickOn(page, el("Run button in \"ANOVA\" dialog")));
    await session.step(34, "Then \"ANOVA\" dialog should be hidden", () => shouldBe(page, el("\"ANOVA\" dialog"), "hidden"));
    await session.step(35, "And box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
    await session.step(36, "And description of box plot viewer should contain text \"affects\"", () => shouldContainText(page, el("description of box plot viewer"), "affects"));
    await session.step(37, "And table \"ANOVA result\" should have 1 row", () => tableRows(page, "ANOVA result", 1));
    await session.step(38, "When the video says \"The answer is on top of the box plot, with its p-value\" at the bottom", () => videoSaysAt(page, "The answer is on top of the box plot, with its p-value", "bottom"));
    await session.step(39, "And user hovers over description of box plot viewer", () => hoverOver(page, el("description of box plot viewer")));
    await session.step(40, "And the video says \"Below: F above F-critical, p below alpha, so the groups differ\" at the top", () => videoSaysAt(page, "Below: F above F-critical, p below alpha, so the groups differ", "top"));
    await session.step(41, "And user hovers over second grid viewer", () => hoverOver(page, el("second grid viewer")));
    await session.step(42, "Then second grid viewer should be bound to table \"ANOVA result\"", () => boundTable(page, el("second grid viewer"), "ANOVA result"));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
