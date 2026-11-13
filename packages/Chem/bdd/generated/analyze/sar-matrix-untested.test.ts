/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-untested.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.sar-matrix]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {matricesBuilt} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix offers the compounds the table holds but never assayed", () => {
  const session = feature(test, "features/analyze/sar-matrix-untested.feature", import.meta.url);
  test("A compound never assayed is offered for testing and queued on the Make list", {tag: ["@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens sar-matrix-para dataset", () => openDataset(page, ds("sar-matrix-para")));
    await session.step(16, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
    await session.step(17, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
    await session.step(18, "When user selects \"smiles\" in Molecules input in \"SAR Matrix\" dialog", () => selectIn(page, "smiles", el("Molecules input in \"SAR Matrix\" dialog")));
    await session.step(19, "And user selects \"pIC50\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "pIC50", el("Activity input in \"SAR Matrix\" dialog")));
    await session.step(20, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")));
    await session.step(21, "And user selects \"Higher is better\" in Direction input in \"SAR Matrix\" dialog", () => selectIn(page, "Higher is better", el("Direction input in \"SAR Matrix\" dialog")));
    await session.step(22, "And user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
    await session.step(23, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(24, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 1", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 1));
    await session.step(25, "And the \"rows\" reading of SAR Matrix Viewer viewer should be 4", () => readingIs(page, "rows", el("SAR Matrix Viewer viewer"), 4));
    await session.step(26, "And the \"columns\" reading of SAR Matrix Viewer viewer should be 4", () => readingIs(page, "columns", el("SAR Matrix Viewer viewer"), 4));
    await session.step(27, "And the \"measured\" reading of SAR Matrix Viewer viewer should be 15", () => readingIs(page, "measured", el("SAR Matrix Viewer viewer"), 15));
    await session.step(28, "And the \"predicted\" reading of SAR Matrix Viewer viewer should be 0", () => readingIs(page, "predicted", el("SAR Matrix Viewer viewer"), 0));
    await session.step(29, "And Never assayed card should contain text \"Already in the table, never assayed · 1\"", () => shouldContainText(page, el("Never assayed card"), "Already in the table, never assayed · 1"));
    await session.step(30, "And first never assayed row should contain text \"n=3\"", () => shouldContainText(page, el("first never assayed row"), "n=3"));
    await session.step(31, "And Queue for testing button should have text \"Queue this compound for testing\"", () => shouldHaveText(page, el("Queue for testing button"), "Queue this compound for testing"));
    await session.step(32, "When user clicks on Queue for testing button", () => clickOn(page, el("Queue for testing button")));
    await session.step(33, "Then an info balloon containing \"Added 1 to the Make list (1 total)\" should have been shown", () => infoBalloonText(page, "Added 1 to the Make list (1 total)"));
    await session.step(34, "And the \"make list\" reading of SAR Matrix Viewer viewer should be 1", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 1));
    await session.step(35, "When user clicks on Make list tab", () => clickOn(page, el("Make list tab")));
    await session.step(36, "Then Make list title should have text \"1 compound\"", () => shouldHaveText(page, el("Make list title"), "1 compound"));
    await session.step(37, "And the \"text of cell 1 of Structure\" reading of Make list grid viewer should be \"Oc1ccc(Cl)cc1\"", () => readingReads(page, "text of cell 1 of Structure", el("Make list grid viewer"), "Oc1ccc(Cl)cc1"));
    await session.step(38, "And the \"text of cell 1 of Status\" reading of Make list grid viewer should be \"Untested\"", () => readingReads(page, "text of cell 1 of Status", el("Make list grid viewer"), "Untested"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
});
