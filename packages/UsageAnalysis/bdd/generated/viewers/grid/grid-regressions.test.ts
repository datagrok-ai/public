/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, public#2348, public#1130, public#1159, public#2718, public#2663, public#629]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, clearSelection, rowsRangeAllSelected, selectNoRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openDatasetRowsAs, openProject, ownProjectGone, saveAllAsProject, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, clickArea, closeContextMenu, doubleClickArea, dragAreaBy, hoverArea, noErrors, noViewerError, pickFromAreaContextMenu, propertyShouldBe, readingHigher, readingLower, readingReads, rightClickArea, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingExcludes, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid regressions from closed issues", () => {
  const session = feature(test, "features/viewers/grid/grid-regressions.feature", import.meta.url);
  test("Sync New Columns decides whether a new calculated column is drawn (public#2348)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#2348"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(22, "When user sets \"Sync New Columns\" property of grid to \"false\"", () => setProperty(page, "Sync New Columns", el("grid"), "false"));
    await session.step(23, "And user adds a calculated column \"H2\" with formula \"${HEIGHT} * 2\"", () => addCalculated(page, "H2", "${HEIGHT} * 2"));
    await session.step(24, "Then the table should have a column \"H2\"", () => hasColumn(page, "H2"));
    await session.step(25, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(26, "When user sets \"Sync New Columns\" property of grid to \"true\"", () => setProperty(page, "Sync New Columns", el("grid"), "true"));
    await session.step(27, "And user adds a calculated column \"H3\" with formula \"${HEIGHT} * 3\"", () => addCalculated(page, "H3", "${HEIGHT} * 3"));
    await session.step(28, "Then the table should have a column \"H3\"", () => hasColumn(page, "H3"));
    await session.step(29, "And the \"column order\" reading of grid should include the text \"H3\"", () => readingIncludes(page, "column order", el("grid"), "H3"));
    await session.step(30, "And the \"column order\" reading of grid should not include the text \"H2\"", () => readingExcludes(page, "column order", el("grid"), "H2"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column narrowed below its name still shows its header tooltip, and the header menu hides it (public#1130)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#1130"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(36, "When user drags the \"column resizer RACE\" area of grid by 60 pixels to the left", () => dragAreaBy(page, "column resizer RACE", el("grid"), 60, "left"));
    await session.step(37, "Then the \"column width of RACE\" reading of grid should be lower than before", () => readingLower(page, "column width of RACE", el("grid")));
    await session.step(38, "When user hovers over the \"header RACE\" area of grid", () => hoverArea(page, "header RACE", el("grid")));
    await session.step(39, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(40, "And tooltip should contain text \"RACE\"", () => shouldContainText(page, el("tooltip"), "RACE"));
    await session.step(41, "When user right-clicks on the \"header RACE\" area of grid", () => rightClickArea(page, "header RACE", el("grid")));
    await session.step(42, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(43, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(44, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A pinned column sorts by its header, descending, ascending and off, and takes a resize (public#1159)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#1159"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(49, "When user picks \"Pin > Pin Column\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header AGE", el("grid")));
    await session.step(50, "Then \"Frozen Columns\" property of grid should be \"2\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "2"));
    await session.step(51, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
    await session.step(52, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
    await session.step(53, "And the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
    await session.step(54, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
    await session.step(55, "Then the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
    await session.step(56, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
    await session.step(57, "Then the \"sort column\" reading of grid should be \"\"", () => readingReads(page, "sort column", el("grid"), ""));
    await session.step(58, "When user drags the \"column resizer AGE\" area of grid by 40 pixels to the right", () => dragAreaBy(page, "column resizer AGE", el("grid"), 40, "right"));
    await session.step(59, "Then the \"column width of AGE\" reading of grid should be higher than before", () => readingHigher(page, "column width of AGE", el("grid")));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
  test("A view with a sorted grid of a second table opens again from its project (public#2718)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#2718"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(65, "Given the user's own project \"bdd grid two tables {time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("bdd grid two tables {time}")));
    await session.step(66, "And user opens demog-1000 dataset keeping the first 100 rows as \"demog-100\"", () => openDatasetRowsAs(page, ds("demog-1000"), 100, "demog-100"));
    await session.step(67, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(68, "And user adds a grid viewer", () => addViewer(page, "grid"));
    await session.step(69, "When user sets \"table\" property of second grid viewer to \"demog-100\"", () => setProperty(page, "table", el("second grid viewer"), "demog-100"));
    await session.step(70, "Then second grid viewer should be bound to table \"demog-100\"", () => boundTable(page, el("second grid viewer"), "demog-100"));
    await session.step(71, "When user double-clicks on the \"header AGE\" area of second grid viewer", () => doubleClickArea(page, "header AGE", el("second grid viewer")));
    await session.step(72, "Then the \"sort column\" reading of second grid viewer should be \"AGE\"", () => readingReads(page, "sort column", el("second grid viewer"), "AGE"));
    await session.step(73, "When user saves all open table views as project \"bdd grid two tables {time}\"", () => saveAllAsProject(page, session.text("bdd grid two tables {time}")));
    await session.step(74, "And user closes all views", () => closeAllViews(page));
    await session.step(75, "And user opens the \"bdd grid two tables {time}\" project", () => openProject(page, session.text("bdd grid two tables {time}")));
    await session.step(76, "Then second grid viewer should be bound to table \"demog-100\"", () => boundTable(page, el("second grid viewer"), "demog-100"));
    await session.step(77, "And the \"sort column\" reading of second grid viewer should be \"AGE\"", () => readingReads(page, "sort column", el("second grid viewer"), "AGE"));
    await session.step(78, "And no viewer of the current view should report an error", () => noViewerError(page));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
  test("Control+Shift+Home and Control+Shift+End select from the current row to the edge (public#2663)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#2663"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(84, "Given user selects no rows", () => selectNoRows(page));
    await session.step(85, "When user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
    await session.step(86, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(87, "When user presses Control+Shift+Home", () => pressKey(page, "Control+Shift+Home"));
    await session.step(88, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
    await session.step(89, "And rows 1 to 5 should all be selected", () => rowsRangeAllSelected(page, 1, 5));
    await session.step(90, "When user clears the row selection", () => clearSelection(page));
    await session.step(91, "And user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
    await session.step(92, "And user presses Control+Shift+End", () => pressKey(page, "Control+Shift+End"));
    await session.step(93, "Then 996 rows should be selected", () => selectedRowCount(page, 996));
    await session.step(94, "And rows 5 to 1000 should all be selected", () => rowsRangeAllSelected(page, 5, 1000));
    await session.step(95, "And no errors should have been logged", () => noErrors(page));
  });
  test("Enter and a typed digit open the editor on the current cell (public#629)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#629"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(100, "When user clicks on the \"cell 4 of AGE\" area of grid", () => clickArea(page, "cell 4 of AGE", el("grid")));
    await session.step(101, "Then cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
    await session.step(102, "When user presses Enter", () => pressKey(page, "Enter"));
    await session.step(103, "Then cell editor should be visible", () => shouldBe(page, el("cell editor"), "visible"));
    await session.step(104, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(105, "Then cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
    await session.step(106, "And the value of \"AGE\" column in row 4 should be \"45\"", () => valueInRow(page, "AGE", 4, "45"));
    await session.step(107, "When user presses 7", () => pressKey(page, "7"));
    await session.step(108, "Then cell editor should be visible", () => shouldBe(page, el("cell editor"), "visible"));
    await session.step(109, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(110, "Then cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
    await session.step(111, "And the value of \"AGE\" column in row 4 should be \"45\"", () => valueInRow(page, "AGE", 4, "45"));
    await session.step(112, "When user presses F2", () => pressKey(page, "F2"));
    await session.step(113, "Then \"AGE\" dialog should be visible", () => shouldBe(page, el("\"AGE\" dialog"), "visible"));
    await session.step(114, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
    await session.step(115, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(116, "Then \"AGE\" dialog should be absent", () => shouldBe(page, el("\"AGE\" dialog"), "absent"));
    await session.step(117, "And no errors should have been logged", () => noErrors(page));
  });
});
