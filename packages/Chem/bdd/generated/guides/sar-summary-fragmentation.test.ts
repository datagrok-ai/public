/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-summary-fragmentation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {callWith} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Summary without fragment columns", () => {
  const session = feature(test, "features/guides/sar-summary-fragmentation.feature", import.meta.url);
  test("no component cards under fragmentation", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And user opens sar-matrix-demo dataset", () => openDataset(page, ds("sar-matrix-demo")));
    await session.step(15, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(16, "When user calls \"Chem:SarMatrixAnalysis\" function with:", () => callWith(page, "Chem:SarMatrixAnalysis", [["table","table"],["molecules","column:smiles"],["activity","column:hERG_pIC50"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""]]), [["table","table"],["molecules","column:smiles"],["activity","column:hERG_pIC50"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""]]);
    await session.step(27, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(28, "When user clicks on \"Summary\" tab", () => clickOn(page, el("\"Summary\" tab")));
    await session.step(29, "Then tab panel should contain text \"Start here\"", () => shouldContainText(page, el("tab panel"), "Start here"));
    await session.step(30, "And tab panel should contain text \"113 L1 · 14 L2\"", () => shouldContainText(page, el("tab panel"), "113 L1 · 14 L2"));
    await session.step(31, "And tab panel should contain text \"open the tab to detect\"", () => shouldContainText(page, el("tab panel"), "open the tab to detect"));
    await session.step(32, "When user clicks on \"SAR Transfer\" tab", () => clickOn(page, el("\"SAR Transfer\" tab")));
    await session.step(33, "Then tab panel should contain text \"transfer\"", () => shouldContainText(page, el("tab panel"), "transfer"));
    await session.step(34, "When user clicks on \"Summary\" tab", () => clickOn(page, el("\"Summary\" tab")));
    await session.step(35, "Then tab panel should not contain text \"open the tab to detect\"", () => shouldNotContainText(page, el("tab panel"), "open the tab to detect"));
    await session.step(36, "When user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(37, "Then tab panel should contain text \"R-groups — within-series ranking\"", () => shouldContainText(page, el("tab panel"), "R-groups — within-series ranking"));
    await session.step(38, "And tab panel should contain text \"Cores are not comparable across series here\"", () => shouldContainText(page, el("tab panel"), "Cores are not comparable across series here"));
    await session.step(39, "And tab panel should not contain text \"offsets from the additive fit\"", () => shouldNotContainText(page, el("tab panel"), "offsets from the additive fit"));
    await session.step(40, "And tab panel should not contain text \"moves pIC50 most\"", () => shouldNotContainText(page, el("tab panel"), "moves pIC50 most"));
    await session.step(41, "When user clicks on \"Method\" summary segment", () => clickOn(page, el("\"Method\" summary segment")));
    await session.step(42, "Then tab panel should contain text \"Fit quality\"", () => shouldContainText(page, el("tab panel"), "Fit quality"));
    await session.step(43, "And tab panel should contain text \"Two different R²\"", () => shouldContainText(page, el("tab panel"), "Two different R²"));
    await session.step(44, "And tab panel should not contain text \"Additive fit over\"", () => shouldNotContainText(page, el("tab panel"), "Additive fit over"));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
