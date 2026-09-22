/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-linked-tables.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addCardFor, linkTables} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, filterIsExactlyCategory, filterPasses, noneOfFiltered, openEmptyFilterPanel, selectFirstRows, selectedRowCount, tableFilterCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel over linked tables", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-linked-tables.feature", import.meta.url);
  test("Filter panel over linked tables", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await session.step(23, "And user opens spgi-linked1 dataset", () => openDataset(page, ds("spgi-linked1")));
    await session.step(24, "And user opens spgi-linked2 dataset", () => openDataset(page, ds("spgi-linked2")));
    await session.step(25, "And the \"SPGI-linked2\" table is linked to the \"SPGI-linked1\" table by \"Sample Name, link column 1, link column 2, link column 3\" to \"Sample Name, link column 1, link column 2, link column 3\" as \"filter to filter\"", () => linkTables(page, "SPGI-linked2", "SPGI-linked1", "Sample Name, link column 1, link column 2, link column 3", "Sample Name, link column 1, link column 2, link column 3", "filter to filter"));
    await session.step(26, "And user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
    await session.step(27, "And user picks \"Data > Link Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Link Tables..."));
    await session.step(28, "And user clicks on \"New Link\" text in \"Link Tables\" dialog", () => clickOn(page, el("\"New Link\" text in \"Link Tables\" dialog")));
    await session.step(29, "And user selects \"selection to filter\" in Link Type input in \"Link Tables\" dialog", () => selectIn(page, "selection to filter", el("Link Type input in \"Link Tables\" dialog")));
    await session.step(30, "And user clicks on LINK button in \"Link Tables\" dialog", () => clickOn(page, el("LINK button in \"Link Tables\" dialog")));
    await session.step(31, "Then \"spgi-100 -> SPGI-linked1\" text in \"Link Tables\" dialog should be visible", () => shouldBe(page, el("\"spgi-100 -> SPGI-linked1\" text in \"Link Tables\" dialog"), "visible"));
    await session.step(32, "When user clicks on CLOSE button in \"Link Tables\" dialog", () => clickOn(page, el("CLOSE button in \"Link Tables\" dialog")));
    await session.step(33, "Then 100 rows of table \"spgi-100\" should pass the filter", () => tableFilterCount(page, 100, "spgi-100"));
    await session.step(34, "And 175 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 175, "SPGI-linked1"));
    await session.step(35, "And 224 rows of table \"SPGI-linked2\" should pass the filter", () => tableFilterCount(page, 224, "SPGI-linked2"));
    await run.scenario("Rows selected in spgi-100 become the filter of SPGI-linked1", async () => {
      await session.step(38, "When user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
      await session.step(39, "And user selects the first 5 rows", () => selectFirstRows(page, 5));
      await session.step(40, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(41, "And 100 rows of table \"spgi-100\" should pass the filter", () => tableFilterCount(page, 100, "spgi-100"));
      await session.step(42, "And 9 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 9, "SPGI-linked1"));
      await session.step(43, "And 224 rows of table \"SPGI-linked2\" should pass the filter", () => tableFilterCount(page, 224, "SPGI-linked2"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card of SPGI-linked2 narrows SPGI-linked1 and nothing comes back", async () => {
      await session.step(47, "When user switches to the \"SPGI-linked2\" table view", () => switchTableView(page, "SPGI-linked2"));
      await session.step(48, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(49, "And user adds a card for \"link column 3\" to the filter panel", () => addCardFor(page, "link column 3"));
      await session.step(50, "And user clicks on the \"category v ii of link column 3\" area of filter panel", () => clickArea(page, "category v ii of link column 3", el("filter panel")));
      await session.step(51, "Then 148 rows should pass the filter", () => filterPasses(page, 148));
      await session.step(52, "And the filter should pass exactly the rows where \"link column 3\" is \"v ii\"", () => filterIsExactlyCategory(page, "link column 3", "v ii"));
      await session.step(53, "And 5 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 5, "SPGI-linked1"));
      await session.step(54, "And 148 rows of table \"SPGI-linked2\" should pass the filter", () => tableFilterCount(page, 148, "SPGI-linked2"));
      await session.step(55, "And 100 rows of table \"spgi-100\" should pass the filter", () => tableFilterCount(page, 100, "spgi-100"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card of SPGI-linked1 composes with what both links bring", async () => {
      await session.step(59, "When user switches to the \"SPGI-linked1\" table view", () => switchTableView(page, "SPGI-linked1"));
      await session.step(60, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(61, "And user adds a card for \"PAMPA Classification\" to the filter panel", () => addCardFor(page, "PAMPA Classification"));
      await session.step(62, "And user clicks on the \"category inconclusive of PAMPA Classification\" area of filter panel", () => clickArea(page, "category inconclusive of PAMPA Classification", el("filter panel")));
      await session.step(63, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(64, "And no rows where \"PAMPA Classification\" is \"> -4.5 cm/s\" should pass the filter", () => noneOfFiltered(page, "PAMPA Classification", "> -4.5 cm/s"));
      await session.step(65, "And no rows where \"PAMPA Classification\" is \"<= -5.3 cm/s\" should pass the filter", () => noneOfFiltered(page, "PAMPA Classification", "<= -5.3 cm/s"));
      await session.step(66, "And 148 rows of table \"SPGI-linked2\" should pass the filter", () => tableFilterCount(page, 148, "SPGI-linked2"));
      await session.step(67, "And 100 rows of table \"spgi-100\" should pass the filter", () => tableFilterCount(page, 100, "spgi-100"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A link changed to selection-to-selection no longer narrows SPGI-linked1", async () => {
      await session.step(71, "When user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
      await session.step(72, "And user clears the row selection", () => clearSelection(page));
      await session.step(73, "Then 51 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 51, "SPGI-linked1"));
      await session.step(74, "When user selects the first 5 rows", () => selectFirstRows(page, 5));
      await session.step(75, "Then 2 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 2, "SPGI-linked1"));
      await session.step(76, "When user picks \"Data > Link Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Link Tables..."));
      await session.step(77, "And user clicks on \"spgi-100 -> SPGI-linked1\" text in \"Link Tables\" dialog", () => clickOn(page, el("\"spgi-100 -> SPGI-linked1\" text in \"Link Tables\" dialog")));
      await session.step(78, "And user selects \"selection to selection\" in Link Type input in \"Link Tables\" dialog", () => selectIn(page, "selection to selection", el("Link Type input in \"Link Tables\" dialog")));
      await session.step(79, "Then Link Type input in \"Link Tables\" dialog should have the value \"selection to selection\"", () => shouldHaveValue(page, el("Link Type input in \"Link Tables\" dialog"), "selection to selection"));
      await session.step(80, "When user clicks on CLOSE button in \"Link Tables\" dialog", () => clickOn(page, el("CLOSE button in \"Link Tables\" dialog")));
      await session.step(81, "Then 51 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 51, "SPGI-linked1"));
      await session.step(82, "And 148 rows of table \"SPGI-linked2\" should pass the filter", () => tableFilterCount(page, 148, "SPGI-linked2"));
      await session.step(83, "When user switches to the \"SPGI-linked1\" table view", () => switchTableView(page, "SPGI-linked1"));
      await session.step(84, "Then 9 rows should be selected", () => selectedRowCount(page, 9));
      await session.step(85, "And 51 rows should pass the filter", () => filterPasses(page, 51));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
