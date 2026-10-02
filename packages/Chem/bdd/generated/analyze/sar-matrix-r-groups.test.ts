/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-r-groups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.sar-matrix]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {matricesBuilt} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingAtLeast, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix from the R-group columns SPGI already has", () => {
  const session = feature(test, "features/analyze/sar-matrix-r-groups.feature", import.meta.url);
  test("SAR Matrix from the R-group columns SPGI already has", {tag: ["@journey", "@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens SPGI-full dataset", () => openDataset(page, ds("SPGI-full")));
    await run.scenario("The dialog offers the R-group columns and states the layout", async () => {
      await session.step(15, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
      await session.step(16, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
      await session.step(17, "And Molecules input in \"SAR Matrix\" dialog should contain text \"Structure\"", () => shouldContainText(page, el("Molecules input in \"SAR Matrix\" dialog"), "Structure"));
      await session.step(18, "When user checks \"Use existing R-groups\" input in \"SAR Matrix\" dialog", () => check(page, el("\"Use existing R-groups\" input in \"SAR Matrix\" dialog")));
      await session.step(19, "Then Core input in \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("Core input in \"SAR Matrix\" dialog"), "visible"));
      await session.step(20, "When user selects \"Core\" in Core input in \"SAR Matrix\" dialog", () => selectIn(page, "Core", el("Core input in \"SAR Matrix\" dialog")));
      await session.step(21, "And user clicks on editor of R-groups input in \"SAR Matrix\" dialog", () => clickOn(page, el("editor of R-groups input in \"SAR Matrix\" dialog")));
      await session.step(22, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(23, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Core\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "Core"));
      await session.step(24, "And the \"text of cell 2 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"R1\"", () => readingReads(page, "text of cell 2 of __name", el("grid viewer in \"Select columns...\" dialog"), "R1"));
      await session.step(25, "And the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"R101\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "R101"));
      await session.step(26, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(27, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(28, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(29, "And user clicks on the \"cell 4 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 4 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(30, "And user clicks on the \"cell 5 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 5 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(31, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(32, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(33, "Then \"Matrix columns\" input in \"SAR Matrix\" dialog should have value \"R2\"", () => shouldHaveValue(page, el("\"Matrix columns\" input in \"SAR Matrix\" dialog"), "R2"));
      await session.step(34, "And \"Matrix columns\" input in \"SAR Matrix\" dialog should offer \"R1, R2, R3, R100, R101\"", () => shouldOffer(page, el("\"Matrix columns\" input in \"SAR Matrix\" dialog"), "R1, R2, R3, R100, R101"));
      await session.step(35, "And \"SAR Matrix\" dialog should contain text \"Rows: Core + R1 + R3 + R100 + R101\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Rows: Core + R1 + R3 + R100 + R101"));
      await session.step(36, "And \"SAR Matrix\" dialog should contain text \"Columns: R2\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Columns: R2"));
      await session.step(37, "When user selects \"Average Mass\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Average Mass", el("Activity input in \"SAR Matrix\" dialog")));
      await session.step(38, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")));
      await session.step(39, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("OK builds one series per core from those columns", async () => {
      await session.step(42, "When user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
      await session.step(43, "Then SAR Matrix Viewer viewer should be visible", () => shouldBe(page, el("SAR Matrix Viewer viewer"), "visible"));
      await session.step(44, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(45, "And the \"source\" reading of SAR Matrix Viewer viewer should be \"R-group columns\"", () => readingReads(page, "source", el("SAR Matrix Viewer viewer"), "R-group columns"));
      await session.step(46, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be at least 5", () => readingAtLeast(page, "matrices", el("SAR Matrix Viewer viewer"), 5));
      await session.step(47, "And the \"predicted with structure\" reading of SAR Matrix Viewer viewer should be at least 1", () => readingAtLeast(page, "predicted with structure", el("SAR Matrix Viewer viewer"), 1));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
