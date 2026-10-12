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
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, downloadContains, downloadCount, downloadThrough, downloadedNotContains, fileDownloaded, shouldBe, uncheck, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickAreaHolding, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A molecule column exported to CSV as SMILES", () => {
  const session = feature(test, "features/general/molecule-in-exported-csv.feature", import.meta.url);
  test("Selected columns exported with Molecules as SMILES", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(14, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(15, "When user clicks on the \"header Id\" area of grid holding Control", () => clickAreaHolding(page, "header Id", el("grid"), "Control"));
    await session.step(16, "And user clicks on the \"header Structure\" area of grid holding Control", () => clickAreaHolding(page, "header Structure", el("grid"), "Control"));
    await session.step(17, "And user clicks on the \"header Chemist\" area of grid holding Control", () => clickAreaHolding(page, "header Chemist", el("grid"), "Control"));
    await session.step(18, "Then columns \"Id, Structure, Chemist\" should be selected", () => columnsSelected(page, "Id, Structure, Chemist"));
    await session.step(21, "Given user watches downloads", () => watchDownloads(page));
    await session.step(22, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(23, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(24, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(25, "When user checks \"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog")));
    await session.step(26, "And user checks \"Selected Columns Only\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Selected Columns Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(27, "And user unchecks \"Selected Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(28, "And user unchecks \"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(29, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(30, "Then a file \"spgi-100.csv\" should have been downloaded", () => fileDownloaded(page, "spgi-100.csv"));
    await session.step(31, "And the downloaded file \"spgi-100.csv\" should contain text \"Id,Structure,Chemist\"", () => downloadContains(page, "spgi-100.csv", "Id,Structure,Chemist"));
    await session.step(32, "And the downloaded file \"spgi-100.csv\" should contain 100 occurrences of \"CAST-\"", () => downloadCount(page, "spgi-100.csv", 100, "CAST-"));
    await session.step(33, "And the downloaded file should not contain \"M  END\"", () => downloadedNotContains(page, "M  END"));
    await session.step(34, "And the downloaded file should not contain \"Last Published Date\"", () => downloadedNotContains(page, "Last Published Date"));
    await session.step(35, "And the downloaded file should not contain \"CAST Idea ID\"", () => downloadedNotContains(page, "CAST Idea ID"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
