/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-r-groups.feature
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
import {check, clickOn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, infoBalloonText, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix from the R-group columns SPGI already has", () => {
  const session = feature(test, "features/analyze/sar-matrix-r-groups.feature", import.meta.url);
  test("SAR Matrix from the R-group columns SPGI already has", {tag: ["@journey", "@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens SPGI-full dataset keeping the first 1500 rows", () => openDatasetRows(page, ds("SPGI-full"), 1500));
    await run.scenario("The dialog says what an R-group set-up lacks and states the layout", async () => {
      await session.step(19, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
      await session.step(20, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
      await session.step(21, "And Molecules input in \"SAR Matrix\" dialog should contain text \"Structure\"", () => shouldContainText(page, el("Molecules input in \"SAR Matrix\" dialog"), "Structure"));
      await session.step(22, "When user checks \"Use existing R-groups\" input in \"SAR Matrix\" dialog", () => check(page, el("\"Use existing R-groups\" input in \"SAR Matrix\" dialog")));
      await session.step(23, "Then \"SAR Matrix\" dialog should contain text \"Pick the core column.\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Pick the core column."));
      await session.step(24, "And OK button in \"SAR Matrix\" dialog should be disabled", () => shouldBe(page, el("OK button in \"SAR Matrix\" dialog"), "disabled"));
      await session.step(25, "When user selects \"Core\" in Core input in \"SAR Matrix\" dialog", () => selectIn(page, "Core", el("Core input in \"SAR Matrix\" dialog")));
      await session.step(26, "Then \"SAR Matrix\" dialog should contain text \"Pick the R-groups.\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Pick the R-groups."));
      await session.step(27, "And OK button in \"SAR Matrix\" dialog should be disabled", () => shouldBe(page, el("OK button in \"SAR Matrix\" dialog"), "disabled"));
      await session.step(28, "When user clicks on editor of R-groups input in \"SAR Matrix\" dialog", () => clickOn(page, el("editor of R-groups input in \"SAR Matrix\" dialog")));
      await session.step(29, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(30, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Core\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "Core"));
      await session.step(31, "And the \"text of cell 2 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"R1\"", () => readingReads(page, "text of cell 2 of __name", el("grid viewer in \"Select columns...\" dialog"), "R1"));
      await session.step(32, "And the \"text of cell 3 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"R2\"", () => readingReads(page, "text of cell 3 of __name", el("grid viewer in \"Select columns...\" dialog"), "R2"));
      await session.step(33, "And the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"R101\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "R101"));
      await session.step(34, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(35, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(36, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(37, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(38, "And user clicks on the \"cell 4 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 4 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(39, "And user clicks on the \"cell 5 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 5 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(40, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(42, "Then an info balloon containing \"Core is the core column, so it was dropped\" should have been shown", () => infoBalloonText(page, "Core is the core column, so it was dropped"));
      await session.step(43, "And \"Matrix columns\" input in \"SAR Matrix\" dialog should have value \"R2\"", () => shouldHaveValue(page, el("\"Matrix columns\" input in \"SAR Matrix\" dialog"), "R2"));
      await session.step(44, "And \"Matrix columns\" input in \"SAR Matrix\" dialog should offer \"R1, R2, R3, R100, R101\"", () => shouldOffer(page, el("\"Matrix columns\" input in \"SAR Matrix\" dialog"), "R1, R2, R3, R100, R101"));
      await session.step(45, "And \"SAR Matrix\" dialog should contain text \"Rows: Core + R1 + R3 + R100 + R101\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Rows: Core + R1 + R3 + R100 + R101"));
      await session.step(46, "And \"SAR Matrix\" dialog should contain text \"Columns: R2\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Columns: R2"));
      await session.step(47, "And OK button in \"SAR Matrix\" dialog should be enabled", () => shouldBe(page, el("OK button in \"SAR Matrix\" dialog"), "enabled"));
      await session.step(48, "When user selects \"Average Mass\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Average Mass", el("Activity input in \"SAR Matrix\" dialog")));
      await session.step(49, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")));
      await session.step(50, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("OK builds one series per core from those columns", async () => {
      await session.step(53, "When user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
      await session.step(54, "Then SAR Matrix Viewer viewer should be visible", () => shouldBe(page, el("SAR Matrix Viewer viewer"), "visible"));
      await session.step(55, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(56, "And the \"source\" reading of SAR Matrix Viewer viewer should be \"R-group columns\"", () => readingReads(page, "source", el("SAR Matrix Viewer viewer"), "R-group columns"));
      await session.step(57, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 6", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 6));
      await session.step(58, "And What to change title should have text \"What to change · 5 components\"", () => shouldHaveText(page, el("What to change title"), "What to change · 5 components"));
      await session.step(59, "When user clicks on more components fold", () => clickOn(page, el("more components fold")));
      await session.step(60, "Then R100 + R101 component row should be visible", () => shouldBe(page, el("R100 + R101 component row"), "visible"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
