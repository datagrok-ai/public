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
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, downloadContains, downloadCount, downloadThrough, downloadedNotContains, fileDownloaded, shouldBe, uncheck, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {columnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickAreaHolding, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A molecule column exported to CSV as SMILES", () => {
  const session = feature(test, "features/general/molecule-in-exported-csv.feature", import.meta.url);
  test("Selected columns exported with Molecules as SMILES", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(16, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(17, "When user clicks on the \"header Id\" area of grid holding Control", () => clickAreaHolding(page, "header Id", el("grid"), "Control"));
    await session.step(18, "And user clicks on the \"header Structure\" area of grid holding Control", () => clickAreaHolding(page, "header Structure", el("grid"), "Control"));
    await session.step(19, "And user clicks on the \"header Chemist\" area of grid holding Control", () => clickAreaHolding(page, "header Chemist", el("grid"), "Control"));
    await session.step(20, "Then columns \"Id, Structure, Chemist\" should be selected", () => columnsSelected(page, "Id, Structure, Chemist"));
    await session.step(23, "Given user watches downloads", () => watchDownloads(page));
    await session.step(24, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(25, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(26, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(27, "When user checks \"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog")));
    await session.step(28, "And user checks \"Selected Columns Only\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Selected Columns Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(29, "And user unchecks \"Selected Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(30, "And user unchecks \"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(31, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(32, "Then a file \"spgi-100.csv\" should have been downloaded", () => fileDownloaded(page, "spgi-100.csv"));
    await session.step(33, "And the downloaded file \"spgi-100.csv\" should contain text \"Id,Structure,Chemist\"", () => downloadContains(page, "spgi-100.csv", "Id,Structure,Chemist"));
    await session.step(34, "And the downloaded file \"spgi-100.csv\" should contain 100 occurrences of \"CAST-\"", () => downloadCount(page, "spgi-100.csv", 100, "CAST-"));
    await session.step(35, "And the downloaded file should not contain \"M  END\"", () => downloadedNotContains(page, "M  END"));
    await session.step(36, "And the downloaded file should not contain \"Last Published Date\"", () => downloadedNotContains(page, "Last Published Date"));
    await session.step(37, "And the downloaded file should not contain \"CAST Idea ID\"", () => downloadedNotContains(page, "CAST Idea ID"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An empty molecule cell is exported as an empty field", {tag: ["@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(16, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(17, "When user clicks on the \"header Id\" area of grid holding Control", () => clickAreaHolding(page, "header Id", el("grid"), "Control"));
    await session.step(18, "And user clicks on the \"header Structure\" area of grid holding Control", () => clickAreaHolding(page, "header Structure", el("grid"), "Control"));
    await session.step(19, "And user clicks on the \"header Chemist\" area of grid holding Control", () => clickAreaHolding(page, "header Chemist", el("grid"), "Control"));
    await session.step(20, "Then columns \"Id, Structure, Chemist\" should be selected", () => columnsSelected(page, "Id, Structure, Chemist"));
    await session.step(43, "Given user opens a table \"mols\" with:", () => openTableOf(page, "mols", [["id","smiles","tag"],["1","CCO","a"],["2","","b"],["3","c1ccccc1","c"]]), [["id","smiles","tag"],["1","CCO","a"],["2","","b"],["3","c1ccccc1","c"]]);
    await session.step(48, "Then \"smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "smiles", "Molecule"));
    await session.step(49, "Given user watches downloads", () => watchDownloads(page));
    await session.step(50, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(51, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(52, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(53, "When user checks \"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog")));
    await session.step(54, "And user unchecks \"Selected Columns Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Columns Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(55, "And user unchecks \"Selected Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(56, "And user unchecks \"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(57, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(58, "Then a file \"mols.csv\" should have been downloaded", () => fileDownloaded(page, "mols.csv"));
    await session.step(59, "And the downloaded file \"mols.csv\" should contain text \"2,,b\"", () => downloadContains(page, "mols.csv", "2,,b"));
    await session.step(60, "And the downloaded file should not contain \"null\"", () => downloadedNotContains(page, "null"));
    await session.step(61, "And the downloaded file \"mols.csv\" should contain text \"1,CCO,a\"", () => downloadContains(page, "mols.csv", "1,CCO,a"));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
});
