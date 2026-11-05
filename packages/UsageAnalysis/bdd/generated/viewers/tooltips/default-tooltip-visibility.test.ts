/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/default-tooltip-visibility.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips, public#2356, public#1599, public#758, public#746]
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
import {clickOn, enterInto, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, colorLinear, colorOff, mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverArea, hoverFirstArea, menuDoesNotList, menuLists, noErrors, openContextMenu, pickFromAreaContextMenu, pickFromContextMenu, pickFromOpenMenu, pointerAway, readingIs, readingReads, readingsDiffer, setProperties, tooltipColumns, tooltipColumnsInOrder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hiding the table's default tooltip from a viewer's context menu", () => {
  const session = feature(test, "features/viewers/tooltips/default-tooltip-visibility.feature", import.meta.url);
  test("Hiding the table's default tooltip from a viewer's context menu", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips", "@realizes:public#2356", "@realizes:public#1599", "@realizes:public#758", "@realizes:public#746"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(36, "Given user is logged in", () => loggedIn(page));
    await session.step(37, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(38, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]);
    await session.step(42, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(43, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(44, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(45, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(46, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(48, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await run.scenario("Before anything is hidden, the grid and both plots show the table's tooltip", async () => {
      await session.step(52, "When user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(53, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(54, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(55, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(56, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(57, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(58, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(59, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(60, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(61, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(62, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(63, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(64, "When user hovers over the \"bin 1\" area of histogram viewer", () => hoverArea(page, "bin 1", el("histogram viewer")));
      await session.step(65, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(66, "And tooltip should contain text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(67, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(68, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(69, "When user hovers over the first \"bar\" area of bar chart viewer", () => hoverFirstArea(page, "bar", el("bar chart viewer")));
      await session.step(70, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(71, "And tooltip should contain text \"count\"", () => shouldContainText(page, el("tooltip"), "count"));
      await session.step(72, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(73, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With a group tooltip set, hovering a histogram bin leaves the grid's colouring and format as they were (public#2356)", async () => {
      await session.step(77, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
      await session.step(78, "And user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(79, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(80, "When user enters \"0.00\" into Custom input in Format AGE dialog", () => enterInto(page, "0.00", el("Custom input in Format AGE dialog")));
      await session.step(81, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(82, "Then Format AGE dialog should be hidden", () => shouldBe(page, el("Format AGE dialog"), "hidden"));
      await session.step(83, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(84, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(85, "When user picks \"Tooltip > Use as Group Tooltip\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Use as Group Tooltip", "view", el("scatter plot viewer")));
      await session.step(86, "And user hovers over the \"bin 8\" area of histogram viewer", () => hoverArea(page, "bin 8", el("histogram viewer")));
      await session.step(87, "Then scatter plot viewer in tooltip should be visible", () => shouldBe(page, el("scatter plot viewer in tooltip"), "visible"));
      await session.step(88, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(89, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(90, "And \"AGE\" column should have tag \"format\" equal to \"0.00\"", () => columnTag(page, "AGE", "format", "0.00"));
      await session.step(91, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(92, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(93, "When user picks \"Tooltip > Remove Group Tooltip\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Remove Group Tooltip", "view", el("scatter plot viewer")));
      await session.step(94, "And user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(95, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tooltip > Hide on the scatter plot hides the row tooltip on the grid and both plots", async () => {
      await session.step(98, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(99, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(100, "When user picks \"Tooltip > Hide\" from the open menu", () => pickFromOpenMenu(page, "Tooltip > Hide"));
      await session.step(101, "And user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(102, "Then the mouse-over row of the table should be 11", () => mouseOverRowIs(page, 11));
      await session.step(103, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(104, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(105, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(106, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(107, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(108, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(109, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(110, "Then the mouse-over row of the table should be 952", () => mouseOverRowIs(page, 952));
      await session.step(111, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(112, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(113, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide switches off the histogram's and the bar chart's own tooltips too", async () => {
      await session.step(116, "When user hovers over the \"bin 1\" area of histogram viewer", () => hoverArea(page, "bin 1", el("histogram viewer")));
      await session.step(117, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(118, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(119, "And user hovers over the first \"bar\" area of bar chart viewer", () => hoverFirstArea(page, "bar", el("bar chart viewer")));
      await session.step(120, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(121, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(122, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Tooltip group offers Show Custom in place of Hide, and Show Custom brings the tooltip back", async () => {
      await session.step(125, "When user opens the context menu of box plot viewer", () => openContextMenu(page, el("box plot viewer")));
      await session.step(126, "Then the open menu should list \"Tooltip > Show Custom\"", () => menuLists(page, "Tooltip > Show Custom"));
      await session.step(127, "And the open menu should not list \"Tooltip > Hide\"", () => menuDoesNotList(page, "Tooltip > Hide"));
      await session.step(128, "When user picks \"Tooltip > Show Custom\" from the open menu", () => pickFromOpenMenu(page, "Tooltip > Show Custom"));
      await session.step(129, "And user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(130, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(131, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(132, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(133, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(134, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(135, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(136, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(137, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(138, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(139, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(140, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(141, "When user hovers over the \"bin 1\" area of histogram viewer", () => hoverArea(page, "bin 1", el("histogram viewer")));
      await session.step(142, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(143, "And tooltip should contain text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(144, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(145, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(146, "When user hovers over the first \"bar\" area of bar chart viewer", () => hoverFirstArea(page, "bar", el("bar chart viewer")));
      await session.step(147, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(148, "And tooltip should contain text \"count\"", () => shouldContainText(page, el("tooltip"), "count"));
      await session.step(149, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(150, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After Design custom tooltip, a column list picked in Edit Tooltip is the tooltip, and Hide hides it (public#758)", async () => {
      await session.step(154, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(155, "And user clicks on \"Design custom tooltip...\" label in \"Edit Tooltip\" dialog", () => clickOn(page, el("\"Design custom tooltip...\" label in \"Edit Tooltip\" dialog")));
      await session.step(156, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(157, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(158, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(159, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(160, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(161, "And user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(162, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(163, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(164, "And user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(165, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
      await session.step(166, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(167, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(168, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(169, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(170, "And the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
      await session.step(171, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(172, "And user picks \"Tooltip > Hide\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Hide", el("scatter plot viewer")));
      await session.step(173, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(174, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(175, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
