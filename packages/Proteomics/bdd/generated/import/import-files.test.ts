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
import {closeAllViews} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {menuFileChooser, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Proteomics > Import opens a MaxQuant proteinGroups.txt as a table of protein groups, directly or through the column mapping", () => {
  const session = feature(test, "features/import/import-files.feature", import.meta.url);
  test("Proteomics > Import opens a MaxQuant proteinGroups.txt as a table of protein groups, directly or through the column mapping", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await run.scenario("MaxQuant... opens the chosen proteinGroups.txt without its contaminant and reverse hits", async () => {
      await session.step(25, "When user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
      await session.step(26, "And user picks \"Proteomics > Import > MaxQuant...\" from the open menu and chooses the \"../files/demo/proteinGroups.txt\" file of the project", () => menuFileChooser(page, "Proteomics > Import > MaxQuant...", "../files/demo/proteinGroups.txt"));
      await session.step(27, "Then the \"proteinGroups\" table view should open with 24 rows", () => tableViewOpened(page, "proteinGroups", 24));
      await session.step(28, "And an info balloon containing \"Imported 24 protein groups\" should have been shown", () => infoBalloonText(page, "Imported 24 protein groups"));
      await session.step(29, "And the value of \"Protein IDs\" column in row 1 should be \"P04406\"", () => valueInRow(page, "Protein IDs", 1, "P04406"));
      await session.step(30, "And \"Protein IDs\" column should have semantic type \"Proteomics-ProteinId\"", () => columnSemType(page, "Protein IDs", "Proteomics-ProteinId"));
      await session.step(31, "And every value of \"Protein IDs\" column should match \"^[A-Z][0-9A-Z]{5}$\"", () => everyValueMatches(page, "Protein IDs", "^[A-Z][0-9A-Z]{5}$"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Generic Matrix... maps the chosen file's columns in a dialog, and OK opens every row", async () => {
      await session.step(35, "When user closes all views", () => closeAllViews(page));
      await session.step(36, "And user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
      await session.step(37, "And user picks \"Proteomics > Import > Generic Matrix...\" from the open menu and chooses the \"../files/demo/proteinGroups.txt\" file of the project", () => menuFileChooser(page, "Proteomics > Import > Generic Matrix...", "../files/demo/proteinGroups.txt"));
      await session.step(38, "Then \"Import Generic Matrix\" dialog should be visible", () => shouldBe(page, el("\"Import Generic Matrix\" dialog"), "visible"));
      await session.step(39, "And editor of \"Protein ID\" input in \"Import Generic Matrix\" dialog should have text \"Protein IDs\"", () => shouldHaveText(page, el("editor of \"Protein ID\" input in \"Import Generic Matrix\" dialog"), "Protein IDs"));
      await session.step(40, "And editor of \"Gene Name (optional)\" input in \"Import Generic Matrix\" dialog should have text \"Gene names\"", () => shouldHaveText(page, el("editor of \"Gene Name (optional)\" input in \"Import Generic Matrix\" dialog"), "Gene names"));
      await session.step(41, "And \"Intensity Columns\" input in \"Import Generic Matrix\" dialog should contain text \"(6) LFQ intensity Sample1\"", () => shouldContainText(page, el("\"Intensity Columns\" input in \"Import Generic Matrix\" dialog"), "(6) LFQ intensity Sample1"));
      await session.step(42, "And \"Log2 Transform\" checkbox in \"Import Generic Matrix\" dialog should be checked", () => shouldBe(page, el("\"Log2 Transform\" checkbox in \"Import Generic Matrix\" dialog"), "checked"));
      await session.step(43, "And \"Import Generic Matrix\" dialog should contain text \"Data appears to be raw intensities\"", () => shouldContainText(page, el("\"Import Generic Matrix\" dialog"), "Data appears to be raw intensities"));
      await session.step(44, "And \"Import Generic Matrix\" dialog should contain text \"Preview: 8 columns, 5 rows\"", () => shouldContainText(page, el("\"Import Generic Matrix\" dialog"), "Preview: 8 columns, 5 rows"));
      await session.step(45, "When user clicks on OK button in \"Import Generic Matrix\" dialog", () => clickOn(page, el("OK button in \"Import Generic Matrix\" dialog")));
      await session.step(46, "Then the \"proteinGroups\" table view should open with 27 rows", () => tableViewOpened(page, "proteinGroups", 27));
      await session.step(47, "And an info balloon containing \"Imported 27 proteins from proteinGroups.txt\" should have been shown", () => infoBalloonText(page, "Imported 27 proteins from proteinGroups.txt"));
      await session.step(48, "And \"Protein IDs\" column should have semantic type \"Proteomics-ProteinId\"", () => columnSemType(page, "Protein IDs", "Proteomics-ProteinId"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
