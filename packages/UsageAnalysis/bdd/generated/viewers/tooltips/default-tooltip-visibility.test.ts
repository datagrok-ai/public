/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/default-tooltip-visibility.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverArea, hoverFirstArea, menuDoesNotList, menuLists, noErrors, openContextMenu, pickFromOpenMenu, pointerAway, readingIs, setProperties, tooltipColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hiding the table's default tooltip from a viewer's context menu", () => {
  const session = feature(test, "features/viewers/tooltips/default-tooltip-visibility.feature", import.meta.url);
  test("Hiding the table's default tooltip from a viewer's context menu", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(35, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]);
    await session.step(39, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(40, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(41, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(42, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(43, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(45, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await run.scenario("Before anything is hidden, the grid and both plots show the table's tooltip", async () => {
      await session.step(49, "When user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(50, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(51, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(52, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(53, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(54, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(55, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(56, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(57, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(58, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(59, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(60, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(61, "When user hovers over the \"bin 1\" area of histogram viewer", () => hoverArea(page, "bin 1", el("histogram viewer")));
      await session.step(62, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(63, "And tooltip should contain text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(64, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(65, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(66, "When user hovers over the first \"bar\" area of bar chart viewer", () => hoverFirstArea(page, "bar", el("bar chart viewer")));
      await session.step(67, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(68, "And tooltip should contain text \"count\"", () => shouldContainText(page, el("tooltip"), "count"));
      await session.step(69, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tooltip > Hide on the scatter plot hides the row tooltip on the grid and both plots", async () => {
      await session.step(73, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(74, "Then the open menu should list \"Tooltip > Hide\"", () => menuLists(page, "Tooltip > Hide"));
      await session.step(75, "When user picks \"Tooltip > Hide\" from the open menu", () => pickFromOpenMenu(page, "Tooltip > Hide"));
      await session.step(76, "And user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(77, "Then the mouse-over row of the table should be 11", () => mouseOverRowIs(page, 11));
      await session.step(78, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(79, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(80, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(81, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(82, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(83, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(84, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(85, "Then the mouse-over row of the table should be 952", () => mouseOverRowIs(page, 952));
      await session.step(86, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(87, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(88, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide switches off the histogram's and the bar chart's own tooltips too", async () => {
      await session.step(91, "When user hovers over the \"bin 1\" area of histogram viewer", () => hoverArea(page, "bin 1", el("histogram viewer")));
      await session.step(92, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(93, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(94, "And user hovers over the first \"bar\" area of bar chart viewer", () => hoverFirstArea(page, "bar", el("bar chart viewer")));
      await session.step(95, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(96, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(97, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Tooltip group offers Show Custom in place of Hide, and Show Custom brings the tooltip back", async () => {
      await session.step(100, "When user opens the context menu of box plot viewer", () => openContextMenu(page, el("box plot viewer")));
      await session.step(101, "Then the open menu should list \"Tooltip > Show Custom\"", () => menuLists(page, "Tooltip > Show Custom"));
      await session.step(102, "And the open menu should not list \"Tooltip > Hide\"", () => menuDoesNotList(page, "Tooltip > Hide"));
      await session.step(103, "When user picks \"Tooltip > Show Custom\" from the open menu", () => pickFromOpenMenu(page, "Tooltip > Show Custom"));
      await session.step(104, "And user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(105, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(106, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(107, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(108, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(109, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(110, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(111, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(112, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(113, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(114, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(115, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(116, "When user hovers over the \"bin 1\" area of histogram viewer", () => hoverArea(page, "bin 1", el("histogram viewer")));
      await session.step(117, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(118, "And tooltip should contain text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(119, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(120, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(121, "When user hovers over the first \"bar\" area of bar chart viewer", () => hoverFirstArea(page, "bar", el("bar chart viewer")));
      await session.step(122, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(123, "And tooltip should contain text \"count\"", () => shouldContainText(page, el("tooltip"), "count"));
      await session.step(124, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(125, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
