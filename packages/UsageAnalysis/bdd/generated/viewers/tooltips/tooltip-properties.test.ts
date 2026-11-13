/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/tooltip-properties.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips, public#2297]
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
import {clickOn, isExpanded, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, hoverArea, menuDoesNotList, menuLists, noErrors, openContextMenu, pickFromContextMenu, pickFromOpenMenu, pointerAway, propertiesShouldBe, propertyShouldBe, readingIs, setProperties, setProperty, tooltipColumns, tooltipColumnsInOrder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A viewer's own tooltip properties", () => {
  const session = feature(test, "features/viewers/tooltips/tooltip-properties.feature", import.meta.url);
  test("A viewer's own tooltip properties", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips", "@realizes:public#2297"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(31, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(32, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(33, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(35, "And user sets properties of second scatter plot viewer:", () => setProperties(page, el("second scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(37, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await run.scenario("Show Tooltip inherits from the table, and Row Tooltip is greyed out and empty", async () => {
      await session.step(41, "When user clicks on settings icon of first scatter plot viewer", () => clickOn(page, el("settings icon of first scatter plot viewer")));
      await session.step(42, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(43, "Given \"Tooltip\" category in context panel is expanded", () => isExpanded(page, el("\"Tooltip\" category in context panel")));
      await session.step(44, "Then \"Show Tooltip\" property in context panel should have value \"inherit from table\"", () => shouldHaveValue(page, el("\"Show Tooltip\" property in context panel"), "inherit from table"));
      await session.step(45, "And \"Row Tooltip\" property in context panel should be disabled", () => shouldBe(page, el("\"Row Tooltip\" property in context panel"), "disabled"));
      await session.step(46, "And \"Row Tooltip\" property of scatter plot viewer should be \"\"", () => propertyShouldBe(page, "Row Tooltip", el("scatter plot viewer"), ""));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom tooltip with no columns of its own lists the viewer's own data columns, and one given columns lists all of them (public#2297)", async () => {
      await session.step(51, "When user selects \"show custom tooltip\" in \"Show Tooltip\" property in context panel", () => selectIn(page, "show custom tooltip", el("\"Show Tooltip\" property in context panel")));
      await session.step(52, "Then \"Show Tooltip\" property of scatter plot viewer should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "show custom tooltip"));
      await session.step(53, "And \"Row Tooltip\" property in context panel should be enabled", () => shouldBe(page, el("\"Row Tooltip\" property in context panel"), "enabled"));
      await session.step(54, "And \"Show Tooltip\" property of second scatter plot viewer should be \"inherit from table\"", () => propertyShouldBe(page, "Show Tooltip", el("second scatter plot viewer"), "inherit from table"));
      await session.step(55, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","HEIGHT"],["Y","WEIGHT"],["Data Values","Merge"]]), [["X","HEIGHT"],["Y","WEIGHT"],["Data Values","Merge"]]);
      await session.step(59, "When user sets \"Show Tooltip\" property of box plot viewer to \"show custom tooltip\"", () => setProperty(page, "Show Tooltip", el("box plot viewer"), "show custom tooltip"));
      await session.step(60, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(61, "Then the mouse-over row of the table should be 952", () => mouseOverRowIs(page, 952));
      await session.step(62, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(63, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(64, "And user sets \"Row Tooltip\" property of box plot viewer to \"AGE\\nSEX\"", () => setProperty(page, "Row Tooltip", el("box plot viewer"), "AGE\\nSEX"));
      await session.step(65, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Column Names","Always"],["Row Tooltip","RACE\nSEX"]]), [["Show Column Names","Always"],["Row Tooltip","RACE\nSEX"]]);
      await session.step(68, "Then \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(69, "When user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(70, "Then the tooltip should show columns \"RACE, SEX\"", () => tooltipColumns(page, "RACE, SEX"));
      await session.step(71, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(72, "Then \"Show Tooltip\" property of second scatter plot viewer should be \"inherit from table\"", () => propertyShouldBe(page, "Show Tooltip", el("second scatter plot viewer"), "inherit from table"));
      await session.step(73, "When user hovers over the \"marker of row 11\" area of second scatter plot viewer", () => hoverArea(page, "marker of row 11", el("second scatter plot viewer")));
      await session.step(74, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(75, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(76, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(77, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(78, "Then the tooltip should show columns \"HEIGHT, WEIGHT\"", () => tooltipColumns(page, "HEIGHT, WEIGHT"));
      await session.step(79, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(80, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(81, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(82, "Then the tooltip should show columns \"AGE, SEX\"", () => tooltipColumns(page, "AGE, SEX"));
      await session.step(83, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(84, "And user sets \"Row Tooltip\" property of scatter plot viewer to \"USUBJID\\nAGE\\nSEX\\nRACE\\nDIS_POP\\nDEMOG\\nCONTROL\\nSTARTED\\nSEVERITY\"", () => setProperty(page, "Row Tooltip", el("scatter plot viewer"), "USUBJID\\nAGE\\nSEX\\nRACE\\nDIS_POP\\nDEMOG\\nCONTROL\\nSTARTED\\nSEVERITY"));
      await session.step(85, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(86, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP, DEMOG, CONTROL, STARTED, SEVERITY\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(87, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(88, "And user sets \"Row Tooltip\" property of scatter plot viewer to \"\"", () => setProperty(page, "Row Tooltip", el("scatter plot viewer"), ""));
      await session.step(89, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tooltip > Hide on a viewer with a custom tooltip hides that tooltip alone", async () => {
      await session.step(92, "When user picks \"Tooltip > Hide\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Hide", el("scatter plot viewer")));
      await session.step(93, "Then \"Show Tooltip\" property of scatter plot viewer should be \"do not show\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "do not show"));
      await session.step(94, "When user opens the context menu of scatter plot viewer", () => openContextMenu(page, el("scatter plot viewer")));
      await session.step(95, "Then the open menu should list \"Tooltip > Show Custom\"", () => menuLists(page, "Tooltip > Show Custom"));
      await session.step(96, "And the open menu should not list \"Tooltip > Hide\"", () => menuDoesNotList(page, "Tooltip > Hide"));
      await session.step(97, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(98, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(99, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(100, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(101, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(102, "And user hovers over the \"marker of row 11\" area of second scatter plot viewer", () => hoverArea(page, "marker of row 11", el("second scatter plot viewer")));
      await session.step(103, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(104, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(105, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(106, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(107, "Then the tooltip should show columns \"AGE, SEX\"", () => tooltipColumns(page, "AGE, SEX"));
      await session.step(108, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(109, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tooltip > Show Custom brings the same custom tooltip back", async () => {
      await session.step(112, "When user picks \"Tooltip > Show Custom\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Show Custom", el("scatter plot viewer")));
      await session.step(113, "Then \"Show Tooltip\" property of scatter plot viewer should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "show custom tooltip"));
      await session.step(114, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(115, "Then the tooltip should show columns \"HEIGHT, WEIGHT\"", () => tooltipColumns(page, "HEIGHT, WEIGHT"));
      await session.step(116, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(117, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"do not show\" in the properties hides the custom tooltip too", async () => {
      await session.step(120, "When user clicks on settings icon of first scatter plot viewer", () => clickOn(page, el("settings icon of first scatter plot viewer")));
      await session.step(121, "Given \"Tooltip\" category in context panel is expanded", () => isExpanded(page, el("\"Tooltip\" category in context panel")));
      await session.step(122, "When user selects \"do not show\" in \"Show Tooltip\" property in context panel", () => selectIn(page, "do not show", el("\"Show Tooltip\" property in context panel")));
      await session.step(123, "Then \"Show Tooltip\" property of scatter plot viewer should be \"do not show\"", () => propertyShouldBe(page, "Show Tooltip", el("scatter plot viewer"), "do not show"));
      await session.step(124, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(125, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(126, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(127, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(128, "And user selects \"show custom tooltip\" in \"Show Tooltip\" property in context panel", () => selectIn(page, "show custom tooltip", el("\"Show Tooltip\" property in context panel")));
      await session.step(129, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(130, "Then the tooltip should show columns \"HEIGHT, WEIGHT\"", () => tooltipColumns(page, "HEIGHT, WEIGHT"));
      await session.step(131, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(132, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide on the reference viewer switches off every tooltip, and Show Custom brings them all back", async () => {
      await session.step(135, "When user picks \"Tooltip > Hide\" from the context menu of second scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Hide", el("second scatter plot viewer")));
      await session.step(136, "And user hovers over the \"marker of row 11\" area of second scatter plot viewer", () => hoverArea(page, "marker of row 11", el("second scatter plot viewer")));
      await session.step(137, "Then the \"hovered row\" reading of second scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("second scatter plot viewer"), 11));
      await session.step(138, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(139, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(140, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(141, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(142, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(143, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(144, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(145, "Then the mouse-over row of the table should be 952", () => mouseOverRowIs(page, 952));
      await session.step(146, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(147, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(148, "And user opens the context menu of second scatter plot viewer", () => openContextMenu(page, el("second scatter plot viewer")));
      await session.step(149, "Then the open menu should list \"Tooltip > Show Custom\"", () => menuLists(page, "Tooltip > Show Custom"));
      await session.step(150, "When user picks \"Tooltip > Show Custom\" from the open menu", () => pickFromOpenMenu(page, "Tooltip > Show Custom"));
      await session.step(151, "And user hovers over the \"marker of row 11\" area of second scatter plot viewer", () => hoverArea(page, "marker of row 11", el("second scatter plot viewer")));
      await session.step(152, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"));
      await session.step(153, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(154, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(155, "Then the tooltip should show columns \"HEIGHT, WEIGHT\"", () => tooltipColumns(page, "HEIGHT, WEIGHT"));
      await session.step(156, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(157, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(158, "Then the tooltip should show columns \"AGE, SEX\"", () => tooltipColumns(page, "AGE, SEX"));
      await session.step(159, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(160, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
