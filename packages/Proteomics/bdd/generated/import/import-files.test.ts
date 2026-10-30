/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/import/import-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, everyValueMatches, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {menuFileChooser, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Proteomics > Import opens a MaxQuant proteinGroups.txt as a table of protein groups, directly or through the column mapping", () => {
  const session = feature(test, "features/import/import-files.feature", import.meta.url);
  test("MaxQuant... opens the chosen proteinGroups.txt without its contaminant and reverse hits", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "When user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
    await session.step(22, "And user picks \"Proteomics > Import > MaxQuant...\" from the open menu and chooses the \"../files/demo/proteinGroups.txt\" file of the project", () => menuFileChooser(page, "Proteomics > Import > MaxQuant...", "../files/demo/proteinGroups.txt"));
    await session.step(23, "Then the \"proteinGroups\" table view should open with 24 rows", () => tableViewOpened(page, "proteinGroups", 24));
    await session.step(24, "And an info balloon containing \"Imported 24 protein groups\" should have been shown", () => infoBalloonText(page, "Imported 24 protein groups"));
    await session.step(25, "And the value of \"Protein IDs\" column in row 1 should be \"P04406\"", () => valueInRow(page, "Protein IDs", 1, "P04406"));
    await session.step(26, "And \"Protein IDs\" column should have semantic type \"Proteomics-ProteinId\"", () => columnSemType(page, "Protein IDs", "Proteomics-ProteinId"));
    await session.step(27, "And every value of \"Protein IDs\" column should match \"^[A-Z][0-9A-Z]{5}$\"", () => everyValueMatches(page, "Protein IDs", "^[A-Z][0-9A-Z]{5}$"));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("Generic Matrix... maps the chosen file's columns in a dialog, and OK opens every row", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "When user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
    await session.step(33, "And user picks \"Proteomics > Import > Generic Matrix...\" from the open menu and chooses the \"../files/demo/proteinGroups.txt\" file of the project", () => menuFileChooser(page, "Proteomics > Import > Generic Matrix...", "../files/demo/proteinGroups.txt"));
    await session.step(34, "Then \"Import Generic Matrix\" dialog should be visible", () => shouldBe(page, el("\"Import Generic Matrix\" dialog"), "visible"));
    await session.step(35, "And editor of \"Protein ID\" input in \"Import Generic Matrix\" dialog should have text \"Protein IDs\"", () => shouldHaveText(page, el("editor of \"Protein ID\" input in \"Import Generic Matrix\" dialog"), "Protein IDs"));
    await session.step(36, "And editor of \"Gene Name (optional)\" input in \"Import Generic Matrix\" dialog should have text \"Gene names\"", () => shouldHaveText(page, el("editor of \"Gene Name (optional)\" input in \"Import Generic Matrix\" dialog"), "Gene names"));
    await session.step(37, "And \"Intensity Columns\" input in \"Import Generic Matrix\" dialog should contain text \"(6) LFQ intensity Sample1\"", () => shouldContainText(page, el("\"Intensity Columns\" input in \"Import Generic Matrix\" dialog"), "(6) LFQ intensity Sample1"));
    await session.step(38, "And \"Log2 Transform\" checkbox in \"Import Generic Matrix\" dialog should be checked", () => shouldBe(page, el("\"Log2 Transform\" checkbox in \"Import Generic Matrix\" dialog"), "checked"));
    await session.step(39, "And \"Import Generic Matrix\" dialog should contain text \"Data appears to be raw intensities\"", () => shouldContainText(page, el("\"Import Generic Matrix\" dialog"), "Data appears to be raw intensities"));
    await session.step(40, "And \"Import Generic Matrix\" dialog should contain text \"Preview: 8 columns, 5 rows\"", () => shouldContainText(page, el("\"Import Generic Matrix\" dialog"), "Preview: 8 columns, 5 rows"));
    await session.step(41, "When user clicks on OK button in \"Import Generic Matrix\" dialog", () => clickOn(page, el("OK button in \"Import Generic Matrix\" dialog")));
    await session.step(42, "Then the \"proteinGroups\" table view should open with 27 rows", () => tableViewOpened(page, "proteinGroups", 27));
    await session.step(43, "And an info balloon containing \"Imported 27 proteins from proteinGroups.txt\" should have been shown", () => infoBalloonText(page, "Imported 27 proteins from proteinGroups.txt"));
    await session.step(44, "And \"Protein IDs\" column should have semantic type \"Proteomics-ProteinId\"", () => columnSemType(page, "Protein IDs", "Proteomics-ProteinId"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
