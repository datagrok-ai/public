/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/tooltips-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips, public#1208, public#2356, public#1599, core#826, public#2297, public#373, public#758, public#746]
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
import {clickOn, enterInto, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, setColumnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedAs, colorLinear} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, dialogCloses, openDataset, openTableOf, switchTableView, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, boundTable, hoverArea, noErrors, pickFromAreaContextMenu, pickFromContextMenu, pointerAway, readingReads, readingsDiffer, setProperties, setProperty, tooltipColumns, tooltipColumnsInOrder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {columnListLacks, columnListStartsWith, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tooltip regressions from closed issues", () => {
  const session = feature(test, "features/viewers/tooltips/tooltips-regressions.feature", import.meta.url);
  test("The edited table tooltip is the same in a cloned view and in the original (public#1208)", {tag: ["@viewers", "@realizes:viewers.tooltips", "@realizes:public#1208"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(24, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(26, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
    await session.step(27, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(28, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(29, "And user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(30, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
    await session.step(31, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(32, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(33, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(34, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
    await session.step(35, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(36, "And user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
    await session.step(37, "Then the \"demog-1000 copy\" view should be current", () => viewIsCurrent(page, "demog-1000 copy"));
    await session.step(38, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(39, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
    await session.step(40, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(41, "And user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
    await session.step(42, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(43, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
    await session.step(44, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(45, "Then no errors should have been logged", () => noErrors(page));
  });
  test("With a group tooltip set, hovering a histogram bin leaves the grid's colouring and format as they were (public#2356)", {tag: ["@viewers", "@realizes:viewers.tooltips", "@realizes:public#2356", "@realizes:public#1599"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(50, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
    await session.step(51, "And user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
    await session.step(52, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
    await session.step(53, "When user enters \"0.00\" into Custom input in Format AGE dialog", () => enterInto(page, "0.00", el("Custom input in Format AGE dialog")));
    await session.step(54, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
    await session.step(55, "Then Format AGE dialog should be hidden", () => shouldBe(page, el("Format AGE dialog"), "hidden"));
    await session.step(56, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
    await session.step(57, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
    await session.step(58, "When user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(59, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(60, "And user picks \"Tooltip > Use as Group Tooltip\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Use as Group Tooltip", "view", el("scatter plot viewer")));
    await session.step(61, "And user hovers over the \"bin 8\" area of histogram viewer", () => hoverArea(page, "bin 8", el("histogram viewer")));
    await session.step(62, "Then scatter plot viewer in tooltip should be visible", () => shouldBe(page, el("scatter plot viewer in tooltip"), "visible"));
    await session.step(63, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(64, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
    await session.step(65, "And \"AGE\" column should have tag \"format\" equal to \"0.00\"", () => columnTag(page, "AGE", "format", "0.00"));
    await session.step(66, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
    await session.step(67, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column drawn by the Color renderer is listed in the row tooltip (core#826)", {tag: ["@viewers", "@realizes:viewers.tooltips", "@realizes:core#826"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(73, "Given the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(74, "And user opens a table \"color-test\" with:", () => openTableOf(page, "color-test", [["id","color","label"],["1","#FF0000","red"],["2","#00FF00","green"],["3","#0000FF","blue"]]), [["id","color","label"],["1","#FF0000","red"],["2","#00FF00","green"],["3","#0000FF","blue"]]);
    await session.step(79, "When user sets the semantic type of \"color\" column to \"Color\"", () => setColumnSemType(page, "color", "Color"));
    await session.step(80, "Then \"color\" column should have semantic type \"Color\"", () => columnSemType(page, "color", "Color"));
    await session.step(81, "And the \"cell type of color\" reading of grid should be \"Color\"", () => readingReads(page, "cell type of color", el("grid"), "Color"));
    await session.step(82, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","id"],["Y","id"],["showLabels","Always"]]), [["X","id"],["Y","id"],["showLabels","Always"]]);
    await session.step(86, "And user hovers over the \"marker of row 2\" area of scatter plot viewer", () => hoverArea(page, "marker of row 2", el("scatter plot viewer")));
    await session.step(87, "Then the tooltip should show columns \"id, color, label\"", () => tooltipColumns(page, "id, color, label"));
    await session.step(88, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(89, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A viewer's custom tooltip lists every column it was given, past the default seven (public#2297)", {tag: ["@viewers", "@realizes:viewers.tooltips", "@realizes:public#2297"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(94, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(95, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"],["Show Tooltip","show custom tooltip"],["Row Tooltip","USUBJID\nAGE\nSEX\nRACE\nDIS_POP\nDEMOG\nCONTROL\nSTARTED\nSEVERITY"]]), [["showLabels","Always"],["Show Tooltip","show custom tooltip"],["Row Tooltip","USUBJID\nAGE\nSEX\nRACE\nDIS_POP\nDEMOG\nCONTROL\nSTARTED\nSEVERITY"]]);
    await session.step(99, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(100, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP, DEMOG, CONTROL, STARTED, SEVERITY\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(101, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(102, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Edit Tooltip on a viewer bound to another table lists that table's columns (public#373)", {tag: ["@viewers", "@realizes:viewers.tooltips", "@realizes:public#373"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(107, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(108, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(109, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(110, "When user sets \"table\" property of scatter plot viewer to \"spgi-100\"", () => setProperty(page, "table", el("scatter plot viewer"), "spgi-100"));
    await session.step(111, "Then scatter plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("scatter plot viewer"), "spgi-100"));
    await session.step(112, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
    await session.step(113, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(114, "And the column list of \"Edit Tooltip\" dialog should start with \"Id\"", () => columnListStartsWith(page, el("\"Edit Tooltip\" dialog"), "Id"));
    await session.step(115, "And the column list of \"Edit Tooltip\" dialog should not list \"USUBJID\"", () => columnListLacks(page, el("\"Edit Tooltip\" dialog"), "USUBJID"));
    await session.step(116, "And the column list of \"Edit Tooltip\" dialog should not list \"AGE\"", () => columnListLacks(page, el("\"Edit Tooltip\" dialog"), "AGE"));
    await session.step(117, "When user clicks on CANCEL button in \"Edit Tooltip\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Tooltip\" dialog")));
    await session.step(118, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(119, "And no errors should have been logged", () => noErrors(page));
  });
  test("After Design custom tooltip, a column list picked in Edit Tooltip is the tooltip, and Hide hides it (public#758)", {tag: ["@viewers", "@realizes:viewers.tooltips", "@realizes:public#758", "@realizes:public#746"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(124, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(125, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(127, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
    await session.step(128, "And user clicks on \"Design custom tooltip...\" label in \"Edit Tooltip\" dialog", () => clickOn(page, el("\"Design custom tooltip...\" label in \"Edit Tooltip\" dialog")));
    await session.step(129, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
    await session.step(130, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(131, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(132, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(133, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(134, "And user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
    await session.step(135, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(136, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
    await session.step(137, "And user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(138, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
    await session.step(139, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(140, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(141, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(142, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
    await session.step(143, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(144, "And user picks \"Tooltip > Hide\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Hide", el("scatter plot viewer")));
    await session.step(145, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
    await session.step(146, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(147, "And no errors should have been logged", () => noErrors(page));
  });
});
