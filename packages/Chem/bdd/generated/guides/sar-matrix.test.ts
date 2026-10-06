/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix", () => {
  const session = feature(test, "features/guides/sar-matrix.feature", import.meta.url);
  test("SAR matrix demo", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens sar-matrix-demo dataset", () => openDataset(page, ds("sar-matrix-demo")));
    await session.step(21, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
    await session.step(22, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
    await session.step(23, "And Activity input in \"SAR Matrix\" dialog should contain text \"CYP3A4\"", () => shouldContainText(page, el("Activity input in \"SAR Matrix\" dialog"), "CYP3A4"));
    await session.step(24, "When user selects \"hERG_pIC50\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "hERG_pIC50", el("Activity input in \"SAR Matrix\" dialog")));
    await session.step(25, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")));
    await session.step(26, "And user selects \"Higher is better\" in Direction input in \"SAR Matrix\" dialog", () => selectIn(page, "Higher is better", el("Direction input in \"SAR Matrix\" dialog")));
    await session.step(27, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(28, "When user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
    await session.step(29, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(30, "And \"Summary\" tab should be visible", () => shouldBe(page, el("\"Summary\" tab"), "visible"));
    await session.step(31, "Then tab panel should contain text \"untransformed\"", () => shouldContainText(page, el("tab panel"), "untransformed"));
    await session.step(32, "And tab panel should contain text \"Cores and substituents found by cutting the molecules\"", () => shouldContainText(page, el("tab panel"), "Cores and substituents found by cutting the molecules"));
    await session.step(33, "And tab panel should contain text \"Start here\"", () => shouldContainText(page, el("tab panel"), "Start here"));
    await session.step(34, "When user clicks on \"L1\" tier chip", () => clickOn(page, el("\"L1\" tier chip")));
    await session.step(35, "Then tab panel should contain text \"Start here\"", () => shouldContainText(page, el("tab panel"), "Start here"));
    await session.step(36, "When user clicks on \"All\" tier chip", () => clickOn(page, el("\"All\" tier chip")));
    await session.step(37, "And user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(38, "Then tab panel should contain text \"within-series ranking\"", () => shouldContainText(page, el("tab panel"), "within-series ranking"));
    await session.step(39, "And tab panel should contain text \"Cores are not comparable across series here\"", () => shouldContainText(page, el("tab panel"), "Cores are not comparable across series here"));
    await session.step(40, "And tab panel should contain text \"What buys potency inside these series\"", () => shouldContainText(page, el("tab panel"), "What buys potency inside these series"));
    await session.step(41, "When user clicks on \"Worth making\" summary segment", () => clickOn(page, el("\"Worth making\" summary segment")));
    await session.step(42, "Then tab panel should contain text \"Already in the table, never assayed\"", () => shouldContainText(page, el("tab panel"), "Already in the table, never assayed"));
    await session.step(43, "And tab panel should contain text \"Worth making\"", () => shouldContainText(page, el("tab panel"), "Worth making"));
    await session.step(44, "When user clicks on \"Method\" summary segment", () => clickOn(page, el("\"Method\" summary segment")));
    await session.step(45, "Then tab panel should contain text \"Fit quality\"", () => shouldContainText(page, el("tab panel"), "Fit quality"));
    await session.step(46, "And tab panel should contain text \"Trust gate\"", () => shouldContainText(page, el("tab panel"), "Trust gate"));
    await session.step(47, "When user clicks on \"Overview\" summary segment", () => clickOn(page, el("\"Overview\" summary segment")));
    await session.step(48, "And user clicks on \"SAR Matrix\" tab", () => clickOn(page, el("\"SAR Matrix\" tab")));
    await session.step(49, "Then \"SAR Matrix\" tab should be visible", () => shouldBe(page, el("\"SAR Matrix\" tab"), "visible"));
    await session.step(50, "When user clicks on \"SAR Transfer\" tab", () => clickOn(page, el("\"SAR Transfer\" tab")));
    await session.step(51, "Then tab panel should not contain text \"Detecting SAR transfers\"", () => shouldNotContainText(page, el("tab panel"), "Detecting SAR transfers"));
    await session.step(52, "When user clicks on \"Make list\" tab", () => clickOn(page, el("\"Make list\" tab")));
    await session.step(53, "And user clicks on \"SAR Matrix\" tab", () => clickOn(page, el("\"SAR Matrix\" tab")));
    await session.step(54, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
