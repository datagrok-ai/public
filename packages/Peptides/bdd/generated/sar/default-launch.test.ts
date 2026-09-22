/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/default-launch.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {analysisScaling} from '../../bindings/activity.js';
import {peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, readingIs, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Launch SAR with the dialog defaults", () => {
  const session = feature(test, "features/sar/default-launch.feature", import.meta.url);
  test("Accepting the defaults builds the analysis on all peptides", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user is logged in", () => loggedIn(page));
    await session.step(7, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(8, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(9, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(10, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(11, "And editor of Sequence input in \"Analyze Peptides\" dialog should have text \"AlignedSequence\"", () => shouldHaveText(page, el("editor of Sequence input in \"Analyze Peptides\" dialog"), "AlignedSequence"));
    await session.step(12, "And editor of Activity input in \"Analyze Peptides\" dialog should have text \"IC50\"", () => shouldHaveText(page, el("editor of Activity input in \"Analyze Peptides\" dialog"), "IC50"));
    await session.step(13, "And Scaling input in \"Analyze Peptides\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"Analyze Peptides\" dialog"), "none"));
    await session.step(14, "And \"Generate clusters\" checkbox in \"Analyze Peptides\" dialog should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in \"Analyze Peptides\" dialog"), "checked"));
    await session.step(15, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(16, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(17, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(18, "And the SAR setting \"mclSettings.threshold\" should be \"70\"", () => sarSetting(page, "mclSettings.threshold", "70"));
    await session.step(19, "And the SAR setting \"mclSettings.inflation\" should be \"1.4\"", () => sarSetting(page, "mclSettings.inflation", "1.4"));
    await session.step(20, "And the \"completed threshold\" reading of MCL viewer should be 70", () => readingIs(page, "completed threshold", el("MCL viewer"), 70));
    await session.step(21, "And the \"completed inflation\" reading of MCL viewer should be 1.4", () => readingIs(page, "completed inflation", el("MCL viewer"), 1.4));
    await session.step(22, "And the SAR activity column should use \"none\" scaling", () => analysisScaling(page, "none"));
    await session.step(23, "And Sequence Variability Map viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Variability Map"));
    await session.step(24, "And Most Potent Residues viewer should be added to the open tableview", () => viewerAdded(page, "Most Potent Residues"));
    await session.step(25, "And Logo Summary Table viewer should be added to the open tableview", () => viewerAdded(page, "Logo Summary Table"));
    await session.step(26, "And the \"members total\" reading of Logo Summary Table viewer should be 647", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 647));
    await session.step(27, "And the \"positions\" reading of Sequence Variability Map viewer should be 17", () => readingIs(page, "positions", el("Sequence Variability Map viewer"), 17));
    await session.step(28, "And scatter plot viewer in MCL viewer should be painted", () => painted(page, el("scatter plot viewer in MCL viewer")));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
