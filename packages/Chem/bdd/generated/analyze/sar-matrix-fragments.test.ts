/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-fragments.feature
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
import {matricesBuilt, readingIsTextOf} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentRowValue, deleteSelected, rowCount, selectFirstRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix cut from the molecules of a small table", () => {
  const session = feature(test, "features/analyze/sar-matrix-fragments.feature", import.meta.url);
  test("SAR Matrix cut from the molecules of a small table", {tag: ["@journey", "@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(15, "And user opens sar-matrix-demo dataset keeping the first 600 rows", () => openDatasetRows(page, ds("sar-matrix-demo"), 600));
    await run.scenario("The dialog cuts the molecules and opens on the Summary", async () => {
      await session.step(18, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
      await session.step(19, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
      await session.step(20, "When user selects \"smiles\" in Molecules input in \"SAR Matrix\" dialog", () => selectIn(page, "smiles", el("Molecules input in \"SAR Matrix\" dialog")));
      await session.step(21, "And user selects \"hERG_pIC50\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "hERG_pIC50", el("Activity input in \"SAR Matrix\" dialog")));
      await session.step(22, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")));
      await session.step(23, "And user selects \"Higher is better\" in Direction input in \"SAR Matrix\" dialog", () => selectIn(page, "Higher is better", el("Direction input in \"SAR Matrix\" dialog")));
      await session.step(24, "And user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
      await session.step(25, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(26, "And the \"source\" reading of SAR Matrix Viewer viewer should be \"fragments\"", () => readingReads(page, "source", el("SAR Matrix Viewer viewer"), "fragments"));
      await session.step(27, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 1", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 1));
      await session.step(28, "And the \"measured\" reading of SAR Matrix Viewer viewer should be 6", () => readingIs(page, "measured", el("SAR Matrix Viewer viewer"), 6));
      await session.step(29, "And the \"predicted\" reading of SAR Matrix Viewer viewer should be 4", () => readingIs(page, "predicted", el("SAR Matrix Viewer viewer"), 4));
      await session.step(30, "And the \"tab\" reading of SAR Matrix Viewer viewer should be \"Summary\"", () => readingReads(page, "tab", el("SAR Matrix Viewer viewer"), "Summary"));
      await session.step(31, "And the \"series\" reading of SAR Matrix Viewer viewer should be the text of name of first key series row", () => readingIsTextOf(page, "series", el("SAR Matrix Viewer viewer"), el("name of first key series row")));
      await session.step(32, "And Worth making card should contain text \"No prediction resting on 2 or more measured compounds on both axes beats its series' best measured compound.\"", () => shouldContainText(page, el("Worth making card"), "No prediction resting on 2 or more measured compounds on both axes beats its series' best measured compound."));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Predictions turned off leave Worth making saying so", async () => {
      await session.step(36, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(37, "And user sets \"Predict virtual analogs\" property of SAR Matrix Viewer viewer to \"false\"", () => setProperty(page, "Predict virtual analogs", el("SAR Matrix Viewer viewer"), "false"));
      await session.step(38, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(39, "And the \"predicted\" reading of SAR Matrix Viewer viewer should be 0", () => readingIs(page, "predicted", el("SAR Matrix Viewer viewer"), 0));
      await session.step(40, "And Worth making card should contain text \"Predicted analogs are off — turn on\"", () => shouldContainText(page, el("Worth making card"), "Predicted analogs are off — turn on"));
      await session.step(41, "When user sets \"Predict virtual analogs\" property of SAR Matrix Viewer viewer to \"true\"", () => setProperty(page, "Predict virtual analogs", el("SAR Matrix Viewer viewer"), "true"));
      await session.step(42, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(43, "And the \"predicted\" reading of SAR Matrix Viewer viewer should be 4", () => readingIs(page, "predicted", el("SAR Matrix Viewer viewer"), 4));
      await session.step(44, "And Worth making card should not contain text \"Predicted analogs are off\"", () => shouldNotContainText(page, el("Worth making card"), "Predicted analogs are off"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rows deleted after a build rebuild the matrices, and a cell still opens its own compound", async () => {
      await session.step(48, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(49, "And user selects the first 10 rows", () => selectFirstRows(page, 10));
      await session.step(50, "And user deletes the selected rows", () => deleteSelected(page));
      await session.step(51, "Then the table should have 590 rows", () => rowCount(page, 590));
      await session.step(52, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
      await session.step(53, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 1", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 1));
      await session.step(54, "And the \"measured\" reading of SAR Matrix Viewer viewer should be 6", () => readingIs(page, "measured", el("SAR Matrix Viewer viewer"), 6));
      await session.step(55, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(56, "And user clicks on the \"cell 1 of c0\" area of grid viewer", () => clickArea(page, "cell 1 of c0", el("grid viewer")));
      await session.step(57, "Then \"CMPD_CHEMBLID\" of the current row should be \"CHEMBL1934525\"", () => currentRowValue(page, "CMPD_CHEMBLID", "CHEMBL1934525"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
