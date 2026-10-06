/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/io/import-export.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.import-export-formats]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, downloadContains, downloadCount, fileDownloaded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Saving a table as SDF", () => {
  const session = feature(test, "features/io/import-export.feature", import.meta.url);
  test("Saving a table as SDF", {tag: ["@journey", "@realizes:chem.cp.import-export-formats"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("Save as SDF offers the molecule column and downloads a record per row", async () => {
      await session.step(16, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(17, "And user watches downloads", () => watchDownloads(page));
      await session.step(18, "When user clicks on \"arrow-to-bottom\" icon in toolbar", () => clickOn(page, el("\"arrow-to-bottom\" icon in toolbar")));
      await session.step(19, "And user clicks on \"As SDF...\" item", () => clickOn(page, el("\"As SDF...\" item")));
      await session.step(20, "Then \"Save as SDF\" dialog should be visible", () => shouldBe(page, el("\"Save as SDF\" dialog"), "visible"));
      await session.step(21, "And Molecules input in \"Save as SDF\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Save as SDF\" dialog"), "canonical_smiles"));
      await session.step(22, "And Notation input in \"Save as SDF\" dialog should have value \"\"", () => shouldHaveValue(page, el("Notation input in \"Save as SDF\" dialog"), ""));
      await session.step(23, "And \"Visible Columns Only\" input in \"Save as SDF\" dialog should be checked", () => shouldBe(page, el("\"Visible Columns Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(24, "And \"Selected Columns Only\" input in \"Save as SDF\" dialog should not be checked", () => shouldNotBe(page, el("\"Selected Columns Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(25, "And \"Filtered Rows Only\" input in \"Save as SDF\" dialog should not be checked", () => shouldNotBe(page, el("\"Filtered Rows Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(26, "And \"Selected Rows Only\" input in \"Save as SDF\" dialog should not be checked", () => shouldNotBe(page, el("\"Selected Rows Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(27, "When user clicks on OK button in \"Save as SDF\" dialog", () => clickOn(page, el("OK button in \"Save as SDF\" dialog")));
      await session.step(28, "Then the \"Save as SDF\" dialog should close", () => dialogCloses(page, "Save as SDF"));
      await session.step(29, "And a file \"smiles.sdf\" should have been downloaded", () => fileDownloaded(page, "smiles.sdf"));
      await session.step(30, "And the downloaded file \"smiles.sdf\" should contain text \"M  END\"", () => downloadContains(page, "smiles.sdf", "M  END"));
      await session.step(31, "And the downloaded file \"smiles.sdf\" should contain text \"V2000\"", () => downloadContains(page, "smiles.sdf", "V2000"));
      await session.step(32, "And the downloaded file \"smiles.sdf\" should contain text \"$$$$\"", () => downloadContains(page, "smiles.sdf", "$$$$"));
      await session.step(33, "And the downloaded file \"smiles.sdf\" should contain text \">  <molregno>\"", () => downloadContains(page, "smiles.sdf", ">  <molregno>"));
      await session.step(34, "And the downloaded file \"smiles.sdf\" should contain 1000 occurrences of \"$$$$\"", () => downloadCount(page, "smiles.sdf", 1000, "$$$$"));
      await session.step(35, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
      await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save as SDF writes V3000 molblocks when that notation is chosen", async () => {
      await session.step(40, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(41, "And user watches downloads", () => watchDownloads(page));
      await session.step(42, "When user clicks on \"arrow-to-bottom\" icon in toolbar", () => clickOn(page, el("\"arrow-to-bottom\" icon in toolbar")));
      await session.step(43, "And user clicks on \"As SDF...\" item", () => clickOn(page, el("\"As SDF...\" item")));
      await session.step(44, "And user selects \"v3Kmolblock\" in Notation input in \"Save as SDF\" dialog", () => selectIn(page, "v3Kmolblock", el("Notation input in \"Save as SDF\" dialog")));
      await session.step(45, "And user clicks on OK button in \"Save as SDF\" dialog", () => clickOn(page, el("OK button in \"Save as SDF\" dialog")));
      await session.step(46, "Then the \"Save as SDF\" dialog should close", () => dialogCloses(page, "Save as SDF"));
      await session.step(47, "And a file \"smiles.sdf\" should have been downloaded", () => fileDownloaded(page, "smiles.sdf"));
      await session.step(48, "And the downloaded file \"smiles.sdf\" should contain text \"V3000\"", () => downloadContains(page, "smiles.sdf", "V3000"));
      await session.step(49, "And the downloaded file \"smiles.sdf\" should contain text \"M  V30 BEGIN ATOM\"", () => downloadContains(page, "smiles.sdf", "M  V30 BEGIN ATOM"));
      await session.step(50, "And the downloaded file \"smiles.sdf\" should contain 1000 occurrences of \"$$$$\"", () => downloadCount(page, "smiles.sdf", 1000, "$$$$"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save as SDF writes only the filtered rows when asked to", async () => {
      await session.step(55, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(56, "And user watches downloads", () => watchDownloads(page));
      await session.step(57, "When user filters rows where \"NumAromaticRings\" is between 2 and 2", () => filterBetween(page, "NumAromaticRings", 2, 2));
      await session.step(58, "Then 254 rows should pass the filter", () => filterPasses(page, 254));
      await session.step(59, "When user clicks on \"arrow-to-bottom\" icon in toolbar", () => clickOn(page, el("\"arrow-to-bottom\" icon in toolbar")));
      await session.step(60, "And user clicks on \"As SDF...\" item", () => clickOn(page, el("\"As SDF...\" item")));
      await session.step(61, "And user checks \"Filtered Rows Only\" input in \"Save as SDF\" dialog", () => check(page, el("\"Filtered Rows Only\" input in \"Save as SDF\" dialog")));
      await session.step(62, "And user clicks on OK button in \"Save as SDF\" dialog", () => clickOn(page, el("OK button in \"Save as SDF\" dialog")));
      await session.step(63, "Then the \"Save as SDF\" dialog should close", () => dialogCloses(page, "Save as SDF"));
      await session.step(64, "And a file \"smiles.sdf\" should have been downloaded", () => fileDownloaded(page, "smiles.sdf"));
      await session.step(65, "And the downloaded file \"smiles.sdf\" should contain 254 occurrences of \"$$$$\"", () => downloadCount(page, "smiles.sdf", 254, "$$$$"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
