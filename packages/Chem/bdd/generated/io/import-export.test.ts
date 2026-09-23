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
import {check, clickOn, downloadContains, downloadCount, downloadFewer, fileDownloaded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnCount, columnSemType, columnUnits, everyValueContains, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPassesFewer, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColors, areaNotColor, areaPainted, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Opening chemical file formats and saving a table as SDF", () => {
  const session = feature(test, "features/io/import-export.feature", import.meta.url);
  test("Opening chemical file formats and saving a table as SDF", {tag: ["@journey", "@realizes:chem.cp.import-export-formats"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("An SDF file opens as one row per record, with its fields beside the molecule", async () => {
      await session.step(18, "Given user opens mol1K.sdf dataset", () => openDataset(page, ds("mol1K.sdf")));
      await session.step(19, "Then the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(20, "And the table should have a column \"molecule\"", () => hasColumn(page, "molecule"));
      await session.step(21, "And \"molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecule", "Molecule"));
      await session.step(22, "And \"molecule\" column should have units \"molblock\"", () => columnUnits(page, "molecule", "molblock"));
      await session.step(23, "And \"molecule\" column should have no missing values", () => columnComplete(page, "molecule"));
      await session.step(24, "And every value of \"molecule\" column should contain \"V2000\"", () => everyValueContains(page, "molecule", "V2000"));
      await session.step(25, "And every value of \"molecule\" column should contain \"M  END\"", () => everyValueContains(page, "molecule", "M  END"));
      await session.step(26, "And the table should have a column \"prID\"", () => hasColumn(page, "prID"));
      await session.step(27, "And the table should have a column \"pIC50_HIV_Integrase\"", () => hasColumn(page, "pIC50_HIV_Integrase"));
      await session.step(28, "And the table should have a column \"Activity_Integrase\"", () => hasColumn(page, "Activity_Integrase"));
      await session.step(29, "And the \"cell 1 of molecule\" area of grid should be painted", () => areaPainted(page, "cell 1 of molecule", el("grid")));
      await session.step(30, "And the \"cell 1 of molecule\" area of grid should be painted in at least 2 colors", () => areaColors(page, "cell 1 of molecule", el("grid"), 2));
      await session.step(31, "And the \"cell 2 of molecule\" area of grid should be painted in at least 2 colors", () => areaColors(page, "cell 2 of molecule", el("grid"), 2));
      await session.step(32, "And the \"cell 1 of prID\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "cell 1 of prID", el("grid"), "#FF0000"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
      await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A MOL2 file opens as one row per TRIPOS molecule", async () => {
      await session.step(37, "Given user opens molecules.mol2 dataset", () => openDataset(page, ds("molecules.mol2")));
      await session.step(38, "Then the table should have 3 rows", () => rowCount(page, 3));
      await session.step(39, "And the table should have 1 column", () => columnCount(page, 1));
      await session.step(40, "And the table should have a column \"molecules\"", () => hasColumn(page, "molecules"));
      await session.step(41, "And \"molecules\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecules", "Molecule"));
      await session.step(42, "And \"molecules\" column should have units \"molblock\"", () => columnUnits(page, "molecules", "molblock"));
      await session.step(43, "And \"molecules\" column should have no missing values", () => columnComplete(page, "molecules"));
      await session.step(44, "And every value of \"molecules\" column should contain \"M  END\"", () => everyValueContains(page, "molecules", "M  END"));
      await session.step(45, "And the \"cell 1 of molecules\" area of grid should be painted", () => areaPainted(page, "cell 1 of molecules", el("grid")));
      await session.step(46, "And the \"cell 1 of molecules\" area of grid should be painted in at least 2 colors", () => areaColors(page, "cell 1 of molecules", el("grid"), 2));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save as SDF offers the molecule column and downloads a record per row", async () => {
      await session.step(51, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(52, "And user watches downloads", () => watchDownloads(page));
      await session.step(53, "When user clicks on \"arrow-to-bottom\" icon in toolbar", () => clickOn(page, el("\"arrow-to-bottom\" icon in toolbar")));
      await session.step(54, "And user clicks on \"As SDF...\" item", () => clickOn(page, el("\"As SDF...\" item")));
      await session.step(55, "Then \"Save as SDF\" dialog should be visible", () => shouldBe(page, el("\"Save as SDF\" dialog"), "visible"));
      await session.step(56, "And Molecules input in \"Save as SDF\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Save as SDF\" dialog"), "canonical_smiles"));
      await session.step(57, "And Notation input in \"Save as SDF\" dialog should have value \"\"", () => shouldHaveValue(page, el("Notation input in \"Save as SDF\" dialog"), ""));
      await session.step(58, "And \"Visible Columns Only\" input in \"Save as SDF\" dialog should be checked", () => shouldBe(page, el("\"Visible Columns Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(59, "And \"Selected Columns Only\" input in \"Save as SDF\" dialog should not be checked", () => shouldNotBe(page, el("\"Selected Columns Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(60, "And \"Filtered Rows Only\" input in \"Save as SDF\" dialog should not be checked", () => shouldNotBe(page, el("\"Filtered Rows Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(61, "And \"Selected Rows Only\" input in \"Save as SDF\" dialog should not be checked", () => shouldNotBe(page, el("\"Selected Rows Only\" input in \"Save as SDF\" dialog"), "checked"));
      await session.step(62, "When user clicks on OK button in \"Save as SDF\" dialog", () => clickOn(page, el("OK button in \"Save as SDF\" dialog")));
      await session.step(63, "Then the \"Save as SDF\" dialog should close", () => dialogCloses(page, "Save as SDF"));
      await session.step(64, "And a file \"smiles.sdf\" should have been downloaded", () => fileDownloaded(page, "smiles.sdf"));
      await session.step(65, "And the downloaded file \"smiles.sdf\" should contain text \"M  END\"", () => downloadContains(page, "smiles.sdf", "M  END"));
      await session.step(66, "And the downloaded file \"smiles.sdf\" should contain text \"V2000\"", () => downloadContains(page, "smiles.sdf", "V2000"));
      await session.step(67, "And the downloaded file \"smiles.sdf\" should contain text \"$$$$\"", () => downloadContains(page, "smiles.sdf", "$$$$"));
      await session.step(68, "And the downloaded file \"smiles.sdf\" should contain text \">  <molregno>\"", () => downloadContains(page, "smiles.sdf", ">  <molregno>"));
      await session.step(69, "And the downloaded file \"smiles.sdf\" should contain 1000 occurrences of \"$$$$\"", () => downloadCount(page, "smiles.sdf", 1000, "$$$$"));
      await session.step(70, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
      await session.step(72, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save as SDF writes V3000 molblocks when that notation is chosen", async () => {
      await session.step(75, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(76, "And user watches downloads", () => watchDownloads(page));
      await session.step(77, "When user clicks on \"arrow-to-bottom\" icon in toolbar", () => clickOn(page, el("\"arrow-to-bottom\" icon in toolbar")));
      await session.step(78, "And user clicks on \"As SDF...\" item", () => clickOn(page, el("\"As SDF...\" item")));
      await session.step(79, "And user selects \"v3Kmolblock\" in Notation input in \"Save as SDF\" dialog", () => selectIn(page, "v3Kmolblock", el("Notation input in \"Save as SDF\" dialog")));
      await session.step(80, "And user clicks on OK button in \"Save as SDF\" dialog", () => clickOn(page, el("OK button in \"Save as SDF\" dialog")));
      await session.step(81, "Then the \"Save as SDF\" dialog should close", () => dialogCloses(page, "Save as SDF"));
      await session.step(82, "And a file \"smiles.sdf\" should have been downloaded", () => fileDownloaded(page, "smiles.sdf"));
      await session.step(83, "And the downloaded file \"smiles.sdf\" should contain text \"V3000\"", () => downloadContains(page, "smiles.sdf", "V3000"));
      await session.step(84, "And the downloaded file \"smiles.sdf\" should contain text \"M  V30 BEGIN ATOM\"", () => downloadContains(page, "smiles.sdf", "M  V30 BEGIN ATOM"));
      await session.step(85, "And the downloaded file \"smiles.sdf\" should contain 1000 occurrences of \"$$$$\"", () => downloadCount(page, "smiles.sdf", 1000, "$$$$"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
      await session.step(87, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save as SDF writes only the filtered rows when asked to", async () => {
      await session.step(90, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(91, "And user watches downloads", () => watchDownloads(page));
      await session.step(92, "When user filters rows where \"NumAromaticRings\" is between 2 and 2", () => filterBetween(page, "NumAromaticRings", 2, 2));
      await session.step(93, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(94, "When user clicks on \"arrow-to-bottom\" icon in toolbar", () => clickOn(page, el("\"arrow-to-bottom\" icon in toolbar")));
      await session.step(95, "And user clicks on \"As SDF...\" item", () => clickOn(page, el("\"As SDF...\" item")));
      await session.step(96, "And user checks \"Filtered Rows Only\" input in \"Save as SDF\" dialog", () => check(page, el("\"Filtered Rows Only\" input in \"Save as SDF\" dialog")));
      await session.step(97, "And user clicks on OK button in \"Save as SDF\" dialog", () => clickOn(page, el("OK button in \"Save as SDF\" dialog")));
      await session.step(98, "Then the \"Save as SDF\" dialog should close", () => dialogCloses(page, "Save as SDF"));
      await session.step(99, "And a file \"smiles.sdf\" should have been downloaded", () => fileDownloaded(page, "smiles.sdf"));
      await session.step(100, "And the downloaded file \"smiles.sdf\" should contain fewer than 1000 occurrences of \"$$$$\"", () => downloadFewer(page, "smiles.sdf", 1000, "$$$$"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
