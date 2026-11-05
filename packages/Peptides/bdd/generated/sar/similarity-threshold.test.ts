/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/similarity-threshold.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {onlyOfSelected, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaAtLeastTall, clickArea, noBalloons, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Launch SAR at a raised similarity threshold", () => {
  const session = feature(test, "features/sar/similarity-threshold.feature", import.meta.url);
  test("Similarity threshold 90 on all peptides preserves the analysis and its selection behavior", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(15, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(16, "Then no errors should have been logged", () => noErrors(page));
    await session.step(17, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(20, "Then the table should have 647 rows", () => rowCount(page, 647));
    await session.step(21, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(22, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(23, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(24, "And user enters \"90\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "90", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(25, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(26, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(27, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(28, "And the SAR setting \"mclSettings.threshold\" should be \"90\"", () => sarSetting(page, "mclSettings.threshold", "90"));
    await session.step(29, "And the \"completed threshold\" reading of MCL viewer should be 90", () => readingIs(page, "completed threshold", el("MCL viewer"), 90));
    await session.step(30, "And the \"members total\" reading of Logo Summary Table viewer should be 647", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 647));
    await session.step(31, "And the \"positions\" reading of Sequence Variability Map viewer should be 17", () => readingIs(page, "positions", el("Sequence Variability Map viewer"), 17));
    await session.step(32, "And the \"header 2\" area of grid should be at least 100 pixels tall", () => areaAtLeastTall(page, "header 2", el("grid"), 100));
    await session.step(33, "When user clicks on the \"A at 2\" area of grid", () => clickArea(page, "A at 2", el("grid")));
    await session.step(34, "Then 299 rows should be selected", () => selectedRowCount(page, 299));
    await session.step(35, "And only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
