/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/molecule-in-exported-csv.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, downloadContains, downloadCount, downloadThrough, downloadedContains, downloadedNotContains, fileDownloaded, shouldBe, uncheck, uploadDownloaded, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnCount, columnSemType, columnUnits, columnsExactly, distinctValues, everyValueMatches, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {columnsSelected, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickAreaHolding, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A molecule column exported to CSV as SMILES reads back as molecules", () => {
  const session = feature(test, "features/general/molecule-in-exported-csv.feature", import.meta.url);
  test("Selected columns exported with Molecules as SMILES, and the file opened again", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "When user clicks on the \"header Id\" area of grid holding Control", () => clickAreaHolding(page, "header Id", el("grid"), "Control"));
    await session.step(17, "And user clicks on the \"header Structure\" area of grid holding Control", () => clickAreaHolding(page, "header Structure", el("grid"), "Control"));
    await session.step(18, "And user clicks on the \"header Chemist\" area of grid holding Control", () => clickAreaHolding(page, "header Chemist", el("grid"), "Control"));
    await session.step(19, "Then columns \"Id, Structure, Chemist\" should be selected", () => columnsSelected(page, "Id, Structure, Chemist"));
    await session.step(22, "Given user watches downloads", () => watchDownloads(page));
    await session.step(23, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(24, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(25, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(26, "When user checks \"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog")));
    await session.step(27, "And user checks \"Selected Columns Only\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Selected Columns Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(28, "And user unchecks \"Selected Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(29, "And user unchecks \"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(30, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(31, "Then a file \"spgi-100.csv\" should have been downloaded", () => fileDownloaded(page, "spgi-100.csv"));
    await session.step(32, "And the downloaded file \"spgi-100.csv\" should contain text \"Id,Structure,Chemist\"", () => downloadContains(page, "spgi-100.csv", "Id,Structure,Chemist"));
    await session.step(33, "And the downloaded file \"spgi-100.csv\" should contain 100 occurrences of \"CAST-\"", () => downloadCount(page, "spgi-100.csv", 100, "CAST-"));
    await session.step(34, "And the downloaded file should not contain \"M  END\"", () => downloadedNotContains(page, "M  END"));
    await session.step(35, "And the downloaded file should not contain \"Last Published Date\"", () => downloadedNotContains(page, "Last Published Date"));
    await session.step(36, "And the downloaded file should not contain \"CAST Idea ID\"", () => downloadedNotContains(page, "CAST Idea ID"));
    await session.step(37, "When user closes all views", () => closeAllViews(page));
    await session.step(38, "And user uploads the downloaded file through \"Open local file\" icon inside browse toolbar", () => uploadDownloaded(page, el("\"Open local file\" icon inside browse toolbar")));
    await session.step(39, "Then the table should have 100 rows", () => rowCount(page, 100));
    await session.step(40, "And the table should have 3 columns", () => columnCount(page, 3));
    await session.step(41, "And the table should have the columns \"Id, Structure, Chemist\"", () => columnsExactly(page, "Id, Structure, Chemist"));
    await session.step(42, "And the value of \"Id\" column in row 1 should be \"CAST-634783\"", () => valueInRow(page, "Id", 1, "CAST-634783"));
    await session.step(43, "And \"Structure\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Structure", "Molecule"));
    await session.step(44, "And \"Structure\" column should have units \"smiles\"", () => columnUnits(page, "Structure", "smiles"));
    await session.step(45, "And every value of \"Structure\" column should match \"^\\S+$\"", () => everyValueMatches(page, "Structure", "^\\S+$"));
    await session.step(46, "And \"Structure\" column should have no missing values", () => columnComplete(page, "Structure"));
    await session.step(47, "And \"Structure\" column should have at least 100 distinct values", () => distinctValues(page, "Structure", 100));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
    await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Without the option the same selection carries the MOLBLOCKs", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "When user clicks on the \"header Id\" area of grid holding Control", () => clickAreaHolding(page, "header Id", el("grid"), "Control"));
    await session.step(17, "And user clicks on the \"header Structure\" area of grid holding Control", () => clickAreaHolding(page, "header Structure", el("grid"), "Control"));
    await session.step(18, "And user clicks on the \"header Chemist\" area of grid holding Control", () => clickAreaHolding(page, "header Chemist", el("grid"), "Control"));
    await session.step(19, "Then columns \"Id, Structure, Chemist\" should be selected", () => columnsSelected(page, "Id, Structure, Chemist"));
    await session.step(52, "Given user watches downloads", () => watchDownloads(page));
    await session.step(53, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(54, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(55, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(56, "When user unchecks \"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog")));
    await session.step(57, "And user checks \"Selected Columns Only\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Selected Columns Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(58, "And user unchecks \"Selected Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(59, "And user unchecks \"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(60, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(61, "Then the downloaded file \"spgi-100.csv\" should contain text \"Id,Structure,Chemist\"", () => downloadContains(page, "spgi-100.csv", "Id,Structure,Chemist"));
    await session.step(62, "And the downloaded file should contain \"M  END\"", () => downloadedContains(page, "M  END"));
    await session.step(63, "And the downloaded file should not contain \"Last Published Date\"", () => downloadedNotContains(page, "Last Published Date"));
    await session.step(64, "And the downloaded file should not contain \"CAST Idea ID\"", () => downloadedNotContains(page, "CAST Idea ID"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
    await session.step(66, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
