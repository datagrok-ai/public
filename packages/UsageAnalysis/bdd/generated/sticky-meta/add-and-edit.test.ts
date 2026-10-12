/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sticky-meta/add-and-edit.feature
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
import {check, clickOn, enterInto, isExpanded, pressKey, shouldBe, shouldContainText, shouldHaveValue, typeAtCaret} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowsRangeSelected, selectFirstRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, dialogCloses, openDataset, packageInstalled, stickyMetaSchemaOnServer, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hoverArea, noBalloons, noErrors, pickFromAreaContextMenu, readingReads, tooltipValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sticky metadata on a molecule cell, as sticky columns and on several rows at once", () => {
  const session = feature(test, "features/sticky-meta/add-and-edit.feature", import.meta.url);
  test("Sticky metadata on a molecule cell, as sticky columns and on several rows at once", {tag: ["@journey", "@serial", "@sticky-meta"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(29, "And a Sticky Meta schema \"bdd-sm-cells-{time}\" for molecules, with the entity type \"bdd-sm-cells-type-{time}\", is on the server with the properties:", () => stickyMetaSchemaOnServer(page, session.text("bdd-sm-cells-{time}"), session.text("bdd-sm-cells-type-{time}"), [["rating","int"],["notes","string"],["verified","bool"],["review_date","datetime"]]), [["rating","int"],["notes","string"],["verified","bool"],["review_date","datetime"]]);
    await session.step(34, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(35, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await run.scenario("Metadata written for one cell shows in its tooltip and stays in its dialog (2.1)", async () => {
      await session.step(38, "When user picks \"Sticky meta > Edit for current cell...\" from the context menu of the \"cell 1 of Structure\" area of grid", () => pickFromAreaContextMenu(page, "Sticky meta > Edit for current cell...", "cell 1 of Structure", el("grid")));
      await session.step(39, "Then \"Sticky meta\" dialog should be visible", () => shouldBe(page, el("\"Sticky meta\" dialog"), "visible"));
      await session.step(40, "When user enters \"5\" into Rating input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog", () => enterInto(page, "5", el(session.text("Rating input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog"))));
      await session.step(41, "And user enters \"test note\" into Notes input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog", () => enterInto(page, "test note", el(session.text("Notes input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog"))));
      await session.step(42, "And user checks Verified input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog", () => check(page, el(session.text("Verified input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog"))));
      await session.step(43, "Then Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog should be enabled", () => shouldBe(page, el(session.text("Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog")), "enabled"));
      await session.step(44, "When user clicks on Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog", () => clickOn(page, el(session.text("Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog"))));
      await session.step(45, "Then Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog should be disabled", () => shouldBe(page, el(session.text("Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog")), "disabled"));
      await session.step(46, "When user clicks on CANCEL button in \"Sticky meta\" dialog", () => clickOn(page, el("CANCEL button in \"Sticky meta\" dialog")));
      await session.step(47, "Then the \"Sticky meta\" dialog should close", () => dialogCloses(page, "Sticky meta"));
      await session.step(48, "When user hovers over the \"top right corner of cell 1 of Structure\" area of grid", () => hoverArea(page, "top right corner of cell 1 of Structure", el("grid")));
      await session.step(49, "Then the tooltip should show \"rating\" as \"5\"", () => tooltipValue(page, "rating", "5"));
      await session.step(50, "And the tooltip should show \"notes\" as \"test note\"", () => tooltipValue(page, "notes", "test note"));
      await session.step(51, "And the tooltip should show \"verified\" as \"true\"", () => tooltipValue(page, "verified", "true"));
      await session.step(52, "When user picks \"Sticky meta > Edit for current cell...\" from the context menu of the \"cell 1 of Structure\" area of grid", () => pickFromAreaContextMenu(page, "Sticky meta > Edit for current cell...", "cell 1 of Structure", el("grid")));
      await session.step(53, "Then Rating input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog should have value \"5\"", () => shouldHaveValue(page, el(session.text("Rating input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog")), "5"));
      await session.step(54, "And Notes input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog should have value \"test note\"", () => shouldHaveValue(page, el(session.text("Notes input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog")), "test note"));
      await session.step(55, "And Verified input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog should be checked", () => shouldBe(page, el(session.text("Verified input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog")), "checked"));
      await session.step(56, "When user clicks on CANCEL button in \"Sticky meta\" dialog", () => clickOn(page, el("CANCEL button in \"Sticky meta\" dialog")));
      await session.step(57, "Then the \"Sticky meta\" dialog should close", () => dialogCloses(page, "Sticky meta"));
      await session.step(58, "When user hovers over the \"top right corner of cell 2 of Structure\" area of grid", () => hoverArea(page, "top right corner of cell 2 of Structure", el("grid")));
      await session.step(59, "Then tooltip should contain text \"No sticky meta for this cell\"", () => shouldContainText(page, el("tooltip"), "No sticky meta for this cell"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The column's Sticky meta pane adds the properties as sticky columns (2.2)", async () => {
      await session.step(64, "When user clicks on the \"header Structure\" area of grid", () => clickArea(page, "header Structure", el("grid")));
      await session.step(65, "Then the context panel should show \"Structure\"", () => contextPanelShows(page, "Structure"));
      await session.step(66, "Given \"Sticky meta\" pane in context panel is expanded", () => isExpanded(page, el("\"Sticky meta\" pane in context panel")));
      await session.step(67, "When user clicks on \"Add bdd-sm-cells-{time}'s properties as columns\" button in context panel", () => clickOn(page, el(session.text("\"Add bdd-sm-cells-{time}'s properties as columns\" button in context panel"))));
      await session.step(68, "Then the table should have a column \"rating\"", () => hasColumn(page, "rating"));
      await session.step(69, "And the table should have a column \"notes\"", () => hasColumn(page, "notes"));
      await session.step(70, "And the table should have a column \"verified\"", () => hasColumn(page, "verified"));
      await session.step(71, "And the table should have a column \"review_date\"", () => hasColumn(page, "review_date"));
      await session.step(72, "And the value of \"rating\" column in row 1 should be \"5\"", () => valueInRow(page, "rating", 1, "5"));
      await session.step(73, "And the value of \"notes\" column in row 1 should be \"test note\"", () => valueInRow(page, "notes", 1, "test note"));
      await session.step(74, "And the value of \"verified\" column in row 1 should be \"true\"", () => valueInRow(page, "verified", 1, "true"));
      await session.step(75, "And the value of \"rating\" column in row 2 should be \"\"", () => valueInRow(page, "rating", 2, ""));
      await session.step(76, "And the value of \"notes\" column in row 2 should be \"\"", () => valueInRow(page, "notes", 2, ""));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
      await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A sticky column sorts like any column (2.2)", async () => {
      await session.step(81, "When user picks \"Sticky meta > Edit for current cell...\" from the context menu of the \"cell 2 of Structure\" area of grid", () => pickFromAreaContextMenu(page, "Sticky meta > Edit for current cell...", "cell 2 of Structure", el("grid")));
      await session.step(82, "And user enters \"9\" into Rating input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog", () => enterInto(page, "9", el(session.text("Rating input in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog"))));
      await session.step(83, "And user clicks on Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog", () => clickOn(page, el(session.text("Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog"))));
      await session.step(84, "Then Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog should be disabled", () => shouldBe(page, el(session.text("Save button in \"bdd-sm-cells-{time}\" section in \"Sticky meta\" dialog")), "disabled"));
      await session.step(85, "When user clicks on CANCEL button in \"Sticky meta\" dialog", () => clickOn(page, el("CANCEL button in \"Sticky meta\" dialog")));
      await session.step(86, "Then the \"Sticky meta\" dialog should close", () => dialogCloses(page, "Sticky meta"));
      await session.step(87, "And the value of \"rating\" column in row 2 should be \"9\"", () => valueInRow(page, "rating", 2, "9"));
      await session.step(88, "When user double-clicks on the \"header rating\" area of grid", () => doubleClickArea(page, "header rating", el("grid")));
      await session.step(89, "Then the \"sort column\" reading of grid should be \"rating\"", () => readingReads(page, "sort column", el("grid"), "rating"));
      await session.step(90, "And the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
      await session.step(91, "And the \"row order\" reading of grid should be \"2, 1, 3, 4, 5, 6, 7, 8, 9, 10\"", () => readingReads(page, "row order", el("grid"), "2, 1, 3, 4, 5, 6, 7, 8, 9, 10"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A removed sticky column comes back with its values (2.2)", async () => {
      await session.step(95, "When user picks \"Remove\" from the context menu of the \"header rating\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header rating", el("grid")));
      await session.step(96, "Then the table should not have a column \"rating\"", () => hasNoColumn(page, "rating"));
      await session.step(97, "When user clicks on the \"header Structure\" area of grid", () => clickArea(page, "header Structure", el("grid")));
      await session.step(98, "Then the context panel should show \"Structure\"", () => contextPanelShows(page, "Structure"));
      await session.step(99, "Given \"Sticky meta\" pane in context panel is expanded", () => isExpanded(page, el("\"Sticky meta\" pane in context panel")));
      await session.step(100, "When user clicks on \"Add rating as a column\" button in context panel", () => clickOn(page, el("\"Add rating as a column\" button in context panel")));
      await session.step(101, "Then the table should have a column \"rating\"", () => hasColumn(page, "rating"));
      await session.step(102, "And the value of \"rating\" column in row 1 should be \"5\"", () => valueInRow(page, "rating", 1, "5"));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
      await session.step(104, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Edit for all properties writes the same values into every selected row (2.3)", async () => {
      await session.step(107, "When user selects the first 3 rows", () => selectFirstRows(page, 3));
      await session.step(108, "Then rows 1 to 3 should be selected", () => rowsRangeSelected(page, 1, 3));
      await session.step(109, "When user picks \"Sticky meta > Edit for all properties\" from the context menu of the \"cell 2 of Structure\" area of grid", () => pickFromAreaContextMenu(page, "Sticky meta > Edit for all properties", "cell 2 of Structure", el("grid")));
      await session.step(110, "Then Rows input should have value \"Selected\"", () => shouldHaveValue(page, el("Rows input"), "Selected"));
      await session.step(111, "When user checks \"verified\" property", () => check(page, el("\"verified\" property")));
      await session.step(112, "And user clicks on editor of \"notes\" property", () => clickOn(page, el("editor of \"notes\" property")));
      await session.step(113, "And user types \"batch note\" at the caret", () => typeAtCaret(page, "batch note"));
      await session.step(114, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(115, "And user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
      await session.step(116, "Then the value of \"notes\" column in row 1 should be \"batch note\"", () => valueInRow(page, "notes", 1, "batch note"));
      await session.step(117, "And the value of \"notes\" column in row 2 should be \"batch note\"", () => valueInRow(page, "notes", 2, "batch note"));
      await session.step(118, "And the value of \"notes\" column in row 3 should be \"batch note\"", () => valueInRow(page, "notes", 3, "batch note"));
      await session.step(119, "And the value of \"verified\" column in row 2 should be \"true\"", () => valueInRow(page, "verified", 2, "true"));
      await session.step(120, "And the value of \"verified\" column in row 3 should be \"true\"", () => valueInRow(page, "verified", 3, "true"));
      await session.step(121, "And the value of \"notes\" column in row 4 should be \"\"", () => valueInRow(page, "notes", 4, ""));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
      await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
