/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/undo-redo.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [public#629]
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
import {clickOn, enterInto, hoverOver, pressKey, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnType, columnUnits, columnsExactly, displayedInRow, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {allOfFiltered, colorCodedAs, filterPasses, filterPassesAll, noColorCoding, resetFilter, tableColumns} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, dialogCloses, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, hoverArea, noErrors, pickFromAreaContextMenu, pointerAway, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A column edited from the table view — removed, renamed, color-coded, its properties and type changed — and Ctrl+Z puts it back", () => {
  const session = feature(test, "features/table-view/undo-redo.feature", import.meta.url);
  test("A column edited from the table view — removed, renamed, color-coded, its properties and type changed — and Ctrl+Z puts it back", {tag: ["@journey", "@realizes:public#629"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens a table \"undo\" with:", () => openTableOf(page, "undo", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]);
    await run.scenario("A column removed from its header menu comes back in its place with Ctrl+Z, and Edit > Redo and Edit > Undo do the same", async () => {
      await session.step(31, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(32, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(33, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(34, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(35, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(36, "And the \"column order\" reading of grid should be \"name, city, weight\"", () => readingReads(page, "column order", el("grid"), "name, city, weight"));
      await session.step(37, "And the value of \"city\" column in row 2 should be \"Rome\"", () => valueInRow(page, "city", 2, "Rome"));
      await session.step(38, "When user picks \"Edit > Redo\" from the top menu", () => pickFromTopMenu(page, "Edit > Redo"));
      await session.step(39, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(40, "And the table should have 2 columns", () => columnCount(page, 2));
      await session.step(41, "When user picks \"Edit > Undo\" from the top menu", () => pickFromTopMenu(page, "Edit > Undo"));
      await session.step(42, "Then the \"column order\" reading of grid should be \"name, city, weight\"", () => readingReads(page, "column order", el("grid"), "name, city, weight"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column renamed in Column Properties gets its name back, and Ctrl+Y renames it again", async () => {
      await session.step(46, "When user picks \"Column Properties...\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header city", el("grid")));
      await session.step(47, "And user enters \"town\" into \"New name:\" input in \"city\" dialog", () => enterInto(page, "town", el("\"New name:\" input in \"city\" dialog")));
      await session.step(48, "And user clicks on OK button in \"city\" dialog", () => clickOn(page, el("OK button in \"city\" dialog")));
      await session.step(49, "Then the table should have a column \"town\"", () => hasColumn(page, "town"));
      await session.step(50, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(51, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(52, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(53, "And the table should not have a column \"town\"", () => hasNoColumn(page, "town"));
      await session.step(54, "And the value of \"city\" column in row 2 should be \"Rome\"", () => valueInRow(page, "city", 2, "Rome"));
      await session.step(55, "When user presses Control+Y", () => pressKey(page, "Control+Y"));
      await session.step(56, "Then the table should have a column \"town\"", () => hasColumn(page, "town"));
      await session.step(57, "And the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(58, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(59, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color coding set from the header menu is taken off by Ctrl+Z and put back by redo", async () => {
      await session.step(63, "When user picks \"Color Coding > Linear\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header weight", el("grid")));
      await session.step(64, "Then \"weight\" column should be color-coded linearly", () => colorCodedAs(page, "weight", "linearly"));
      await session.step(65, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(66, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(67, "Then \"weight\" column should have no color coding", () => noColorCoding(page, "weight"));
      await session.step(68, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(69, "Then \"weight\" column should be color-coded linearly", () => colorCodedAs(page, "weight", "linearly"));
      await session.step(70, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(71, "Then \"weight\" column should have no color coding", () => noColorCoding(page, "weight"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Shift+R clears the view down to the grid, Ctrl+Z brings the viewers back and redo clears again", async () => {
      await session.step(75, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(76, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(77, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(78, "And user presses Control+Shift+R", () => pressKey(page, "Control+Shift+R"));
      await session.step(79, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(80, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(81, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(82, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(83, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(84, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(85, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(86, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(87, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("F2 opens the current column's Column Properties, not the cell editor, and a format tag set there changes the grid until Ctrl+Z (public#629)", async () => {
      await session.step(93, "Then the \"weight\" cell of row 1 should be displayed as \"1.50\"", () => displayedInRow(page, "weight", 1, "1.50"));
      await session.step(94, "When user clicks on the \"cell 1 of weight\" area of grid", () => clickArea(page, "cell 1 of weight", el("grid")));
      await session.step(95, "And user presses F2", () => pressKey(page, "F2"));
      await session.step(96, "Then \"weight\" dialog should be visible", () => shouldBe(page, el("\"weight\" dialog"), "visible"));
      await session.step(97, "And cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
      await session.step(98, "And \"New name:\" input in \"weight\" dialog should have value \"weight\"", () => shouldHaveValue(page, el("\"New name:\" input in \"weight\" dialog"), "weight"));
      await session.step(99, "And \"Key 1\" text input in \"weight\" dialog should have value \"format\"", () => shouldHaveValue(page, el("\"Key 1\" text input in \"weight\" dialog"), "format"));
      await session.step(100, "And \"Value 1\" text input in \"weight\" dialog should have value \"0.00\"", () => shouldHaveValue(page, el("\"Value 1\" text input in \"weight\" dialog"), "0.00"));
      await session.step(101, "When user enters \"0.000\" into \"Value 1\" text input in \"weight\" dialog", () => enterInto(page, "0.000", el("\"Value 1\" text input in \"weight\" dialog")));
      await session.step(102, "And user clicks on OK button in \"weight\" dialog", () => clickOn(page, el("OK button in \"weight\" dialog")));
      await session.step(103, "Then the \"weight\" dialog should close", () => dialogCloses(page, "weight"));
      await session.step(104, "And the \"weight\" cell of row 1 should be displayed as \"1.500\"", () => displayedInRow(page, "weight", 1, "1.500"));
      await session.step(105, "And the \"weight\" cell of row 2 should be displayed as \"2.250\"", () => displayedInRow(page, "weight", 2, "2.250"));
      await session.step(106, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(107, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(108, "Then the \"weight\" cell of row 1 should be displayed as \"1.50\"", () => displayedInRow(page, "weight", 1, "1.50"));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Description and units added as two tag rows show in the header tooltip", async () => {
      await session.step(112, "When user picks \"Column Properties...\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header weight", el("grid")));
      await session.step(113, "And user clicks on first \"Add row\" button in \"weight\" dialog", () => clickOn(page, el("first \"Add row\" button in \"weight\" dialog")));
      await session.step(114, "And user enters \"description\" into \"Key 1\" text input in \"weight\" dialog", () => enterInto(page, "description", el("\"Key 1\" text input in \"weight\" dialog")));
      await session.step(115, "And user enters \"Body weight\" into \"Value 1\" text input in \"weight\" dialog", () => enterInto(page, "Body weight", el("\"Value 1\" text input in \"weight\" dialog")));
      await session.step(116, "And user clicks on first \"Add row\" button in \"weight\" dialog", () => clickOn(page, el("first \"Add row\" button in \"weight\" dialog")));
      await session.step(117, "Then \"Key 2\" text input in \"weight\" dialog should have value \"description\"", () => shouldHaveValue(page, el("\"Key 2\" text input in \"weight\" dialog"), "description"));
      await session.step(118, "When user enters \"units\" into \"Key 1\" text input in \"weight\" dialog", () => enterInto(page, "units", el("\"Key 1\" text input in \"weight\" dialog")));
      await session.step(119, "And user enters \"kg\" into \"Value 1\" text input in \"weight\" dialog", () => enterInto(page, "kg", el("\"Value 1\" text input in \"weight\" dialog")));
      await session.step(120, "And user clicks on OK button in \"weight\" dialog", () => clickOn(page, el("OK button in \"weight\" dialog")));
      await session.step(121, "Then the \"weight\" dialog should close", () => dialogCloses(page, "weight"));
      await session.step(122, "And \"weight\" column should have units \"kg\"", () => columnUnits(page, "weight", "kg"));
      await session.step(123, "When user hovers over the \"header weight\" area of grid", () => hoverArea(page, "header weight", el("grid")));
      await session.step(124, "Then tooltip should contain text \"Body weight\"", () => shouldContainText(page, el("tooltip"), "Body weight"));
      await session.step(125, "And tooltip should contain text \"units: kg\"", () => shouldContainText(page, el("tooltip"), "units: kg"));
      await session.step(126, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(127, "And user hovers over the \"header city\" area of grid", () => hoverArea(page, "header city", el("grid")));
      await session.step(128, "Then tooltip should contain text \"city\"", () => shouldContainText(page, el("tooltip"), "city"));
      await session.step(129, "And tooltip should not contain text \"Body weight\"", () => shouldNotContainText(page, el("tooltip"), "Body weight"));
      await session.step(130, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(131, "And user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(132, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(133, "Then \"weight\" column should have units \"\"", () => columnUnits(page, "weight", ""));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column Properties refuses a name another column already has", async () => {
      await session.step(137, "When user picks \"Column Properties...\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header weight", el("grid")));
      await session.step(138, "And user enters \"City\" into \"New name:\" input in \"weight\" dialog", () => enterInto(page, "City", el("\"New name:\" input in \"weight\" dialog")));
      await session.step(139, "Then \"New name:\" input in \"weight\" dialog should be invalid", () => shouldBe(page, el("\"New name:\" input in \"weight\" dialog"), "invalid"));
      await session.step(140, "When user hovers over \"New name:\" input in \"weight\" dialog", () => hoverOver(page, el("\"New name:\" input in \"weight\" dialog")));
      await session.step(141, "Then tooltip should contain text \"Column 'City' already exists.\"", () => shouldContainText(page, el("tooltip"), "Column 'City' already exists."));
      await session.step(142, "When user enters \"mass\" into \"New name:\" input in \"weight\" dialog", () => enterInto(page, "mass", el("\"New name:\" input in \"weight\" dialog")));
      await session.step(143, "Then \"New name:\" input in \"weight\" dialog should be valid", () => shouldBe(page, el("\"New name:\" input in \"weight\" dialog"), "valid"));
      await session.step(144, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(145, "Then \"weight\" dialog should be hidden", () => shouldBe(page, el("\"weight\" dialog"), "hidden"));
      await session.step(146, "And the table should have the columns \"name, city, weight\"", () => columnsExactly(page, "name, city, weight"));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("While a dialog is open Ctrl+Z leaves the table alone, and after it closes Ctrl+Z works", async () => {
      await session.step(150, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(151, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(152, "When user clicks on the \"cell 1 of weight\" area of grid", () => clickArea(page, "cell 1 of weight", el("grid")));
      await session.step(153, "And user presses F2", () => pressKey(page, "F2"));
      await session.step(154, "Then \"weight\" dialog should be visible", () => shouldBe(page, el("\"weight\" dialog"), "visible"));
      await session.step(155, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(156, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(157, "And \"weight\" dialog should be visible", () => shouldBe(page, el("\"weight\" dialog"), "visible"));
      await session.step(158, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(159, "Then \"weight\" dialog should be hidden", () => shouldBe(page, el("\"weight\" dialog"), "hidden"));
      await session.step(160, "And the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(161, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(162, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(163, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Z in another table's view does not undo this table's change, and switching back brings the history", async () => {
      await session.step(167, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(168, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(169, "Given user opens a table \"other\" with:", () => openTableOf(page, "other", [["id","score"],["1","10"],["2","20"]]), [["id","score"],["1","10"],["2","20"]]);
      await session.step(173, "When user clicks on the \"cell 1 of id\" area of grid", () => clickArea(page, "cell 1 of id", el("grid")));
      await session.step(174, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(175, "Then table \"undo\" should have columns \"name, weight\"", () => tableColumns(page, "undo", "name, weight"));
      await session.step(176, "And table \"other\" should have columns \"id, score\"", () => tableColumns(page, "other", "id, score"));
      await session.step(177, "When user switches to the \"undo\" table view", () => switchTableView(page, "undo"));
      await session.step(178, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(179, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(180, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(181, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(182, "And table \"other\" should have columns \"id, score\"", () => tableColumns(page, "other", "id, score"));
      await session.step(183, "When user switches to the \"other\" table view", () => switchTableView(page, "other"));
      await session.step(184, "And user closes the current view", () => closeCurrentView(page));
      await session.step(185, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Change Type... in the header menu shows the current type, and OK converts the column, a value that cannot be converted taking the Not convertible value", async () => {
      await session.step(188, "Given user opens a table \"ct\" with:", () => openTableOf(page, "ct", [["name","code"],["a","10"],["b","20"],["c","n/a"],["d",""]]), [["name","code"],["a","10"],["b","20"],["c","n/a"],["d",""]]);
      await session.step(194, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(195, "When user picks \"Change Type...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Change Type...", "header code", el("grid")));
      await session.step(196, "Then \"Change column type\" dialog should be visible", () => shouldBe(page, el("\"Change column type\" dialog"), "visible"));
      await session.step(197, "And \"Type\" input in \"Change column type\" dialog should have value \"string\"", () => shouldHaveValue(page, el("\"Type\" input in \"Change column type\" dialog"), "string"));
      await session.step(198, "And \"Type\" input in \"Change column type\" dialog should be disabled", () => shouldBe(page, el("\"Type\" input in \"Change column type\" dialog"), "disabled"));
      await session.step(199, "When user selects \"int\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "int", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(200, "And user enters \"-1\" into \"Not convertible\" input in \"Change column type\" dialog", () => enterInto(page, "-1", el("\"Not convertible\" input in \"Change column type\" dialog")));
      await session.step(201, "And user clicks on OK button in \"Change column type\" dialog", () => clickOn(page, el("OK button in \"Change column type\" dialog")));
      await session.step(202, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(203, "And \"code\" column should have type \"int\"", () => columnType(page, "code", "int"));
      await session.step(204, "And the value of \"code\" column in row 2 should be \"20\"", () => valueInRow(page, "code", 2, "20"));
      await session.step(205, "And the value of \"code\" column in row 3 should be \"-1\"", () => valueInRow(page, "code", 3, "-1"));
      await session.step(206, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(207, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(208, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(209, "And the value of \"code\" column in row 3 should be \"n/a\"", () => valueInRow(page, "code", 3, "n/a"));
      await session.step(210, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter icon beside Not convertible keeps only the rows that cannot be converted", async () => {
      await session.step(213, "When user picks \"Change Type...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Change Type...", "header code", el("grid")));
      await session.step(214, "And user selects \"int\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "int", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(215, "And user clicks on second \"Filter\" icon in \"Change column type\" dialog", () => clickOn(page, el("second \"Filter\" icon in \"Change column type\" dialog")));
      await session.step(216, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(217, "And all rows where \"code\" is \"n/a\" should pass the filter", () => allOfFiltered(page, "code", "n/a"));
      await session.step(218, "When user clicks on CANCEL button in \"Change column type\" dialog", () => clickOn(page, el("CANCEL button in \"Change column type\" dialog")));
      await session.step(219, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(220, "When user resets the filter", () => resetFilter(page));
      await session.step(221, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(222, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Change type... in Column Properties opens the same dialog, and its OK converts the column", async () => {
      await session.step(225, "When user picks \"Column Properties...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header code", el("grid")));
      await session.step(226, "Then \"code\" dialog should be visible", () => shouldBe(page, el("\"code\" dialog"), "visible"));
      await session.step(227, "When user clicks on \"Change type...\" label in \"code\" dialog", () => clickOn(page, el("\"Change type...\" label in \"code\" dialog")));
      await session.step(228, "Then \"Change column type\" dialog should be visible", () => shouldBe(page, el("\"Change column type\" dialog"), "visible"));
      await session.step(229, "When user selects \"double\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "double", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(230, "And user clicks on OK button in \"Change column type\" dialog", () => clickOn(page, el("OK button in \"Change column type\" dialog")));
      await session.step(231, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(232, "And \"code\" column should have type \"double\"", () => columnType(page, "code", "double"));
      await session.step(233, "And the value of \"code\" column in row 1 should be \"10\"", () => valueInRow(page, "code", 1, "10"));
      await session.step(234, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(235, "Then \"code\" dialog should be hidden", () => shouldBe(page, el("\"code\" dialog"), "hidden"));
      await session.step(236, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
