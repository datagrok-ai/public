/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/grid-visible-columns-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips, core#441, core#407, core#474]
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
import {check, clickOn, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, shouldNotContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, dialogCloses, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaBy, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, propertyShouldBe, readingHigherThanRemembered, readingReads, rememberReading, setProperties, tooltipColumns, tooltipValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The grid's tooltip and the columns it cannot show", () => {
  const session = feature(test, "features/viewers/tooltips/grid-visible-columns-tooltip.feature", import.meta.url);
  test("The grid's tooltip and the columns it cannot show", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips", "@realizes:core#441", "@realizes:core#407", "@realizes:core#474"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
    await run.scenario("By default the grid shows no row tooltip, even for the columns pushed off it", async () => {
      await session.step(36, "Then \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(37, "And \"Row Tooltip\" property of grid should be \"\"", () => propertyShouldBe(page, "Row Tooltip", el("grid"), ""));
      await session.step(38, "And \"Show Visible Columns In Tooltip\" property of grid should be \"false\"", () => propertyShouldBe(page, "Show Visible Columns In Tooltip", el("grid"), "false"));
      await session.step(39, "When user drags the \"column resizer value\" area of grid by 1000 pixels to the right", () => dragAreaBy(page, "column resizer value", el("grid"), 1000, "right"));
      await session.step(40, "Then grid should not have a \"cell 1 of target\" area", () => hasNoArea(page, el("grid"), "cell 1 of target"));
      await session.step(41, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(42, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(43, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(44, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(45, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("By default a cut-off value shows in full, a value that fits shows nothing, and a Row Tooltip shows only itself (core#441)", async () => {
      await session.step(49, "When user closes all views", () => closeAllViews(page));
      await session.step(50, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
      await session.step(57, "And user drags the \"column resizer source\" area of grid by 80 pixels to the left", () => dragAreaBy(page, "column resizer source", el("grid"), 80, "left"));
      await session.step(58, "And user hovers over the \"cell 1 of source\" area of grid", () => hoverArea(page, "cell 1 of source", el("grid")));
      await session.step(59, "Then tooltip should contain text \"Agricultural 'waste'\"", () => shouldContainText(page, el("tooltip"), "Agricultural 'waste'"));
      await session.step(60, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(61, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(62, "When user hovers over the \"cell 1 of target\" area of grid", () => hoverArea(page, "cell 1 of target", el("grid")));
      await session.step(63, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(64, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(65, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Row Height","80"]]), [["Row Height","80"]]);
      await session.step(67, "And user hovers over the \"cell 1 of source\" area of grid", () => hoverArea(page, "cell 1 of source", el("grid")));
      await session.step(68, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(69, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(70, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(71, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Row Height","20"],["Show Column Names","Always"],["Row Tooltip","target"]]), [["Row Height","20"],["Show Column Names","Always"],["Row Tooltip","target"]]);
      await session.step(75, "And user hovers over the \"cell 1 of source\" area of grid", () => hoverArea(page, "cell 1 of source", el("grid")));
      await session.step(76, "Then the tooltip should show columns \"target\"", () => tooltipColumns(page, "target"));
      await session.step(77, "And tooltip should not contain text \"Agricultural\"", () => shouldNotContainText(page, el("tooltip"), "Agricultural"));
      await session.step(78, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(79, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("View > Edit Tooltip gives the grid the columns picked, unless it has its own (core#474)", async () => {
      await session.step(83, "When user closes all views", () => closeAllViews(page));
      await session.step(84, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
      await session.step(91, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Column Names","Always"]]), [["Show Column Names","Always"]]);
      await session.step(93, "And user picks \"View > Edit Tooltip...\" from the top menu", () => pickFromTopMenu(page, "View > Edit Tooltip..."));
      await session.step(94, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(95, "And \"Use tooltip\" input in \"Edit Tooltip\" dialog should have value \"Table\"", () => shouldHaveValue(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "Table"));
      await session.step(96, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(97, "And user toggles the \"target\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "target", el("\"Edit Tooltip\" dialog")));
      await session.step(98, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(99, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(100, "And \"Show Tooltip\" property of grid should be \"inherit from table\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "inherit from table"));
      await session.step(101, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(102, "Then the tooltip should show columns \"target\"", () => tooltipColumns(page, "target"));
      await session.step(103, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(104, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","show custom tooltip"],["Row Tooltip","source"]]), [["Show Tooltip","show custom tooltip"],["Row Tooltip","source"]]);
      await session.step(107, "And user picks \"View > Edit Tooltip...\" from the top menu", () => pickFromTopMenu(page, "View > Edit Tooltip..."));
      await session.step(108, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(109, "And \"Use tooltip\" input in \"Edit Tooltip\" dialog should have value \"Grid\"", () => shouldHaveValue(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "Grid"));
      await session.step(110, "When user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(111, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(112, "And \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(113, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(114, "Then the tooltip should show columns \"source\"", () => tooltipColumns(page, "source"));
      await session.step(115, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(116, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("View > Edit Tooltip can give a default grid its own Row Tooltip (core#474)", async () => {
      await session.step(120, "When user closes all views", () => closeAllViews(page));
      await session.step(121, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
      await session.step(128, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Column Names","Always"]]), [["Show Column Names","Always"]]);
      await session.step(130, "And user picks \"View > Edit Tooltip...\" from the top menu", () => pickFromTopMenu(page, "View > Edit Tooltip..."));
      await session.step(131, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(132, "When user selects \"Grid\" in \"Use tooltip\" input in \"Edit Tooltip\" dialog", () => selectIn(page, "Grid", el("\"Use tooltip\" input in \"Edit Tooltip\" dialog")));
      await session.step(133, "And user toggles the \"target\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "target", el("\"Edit Tooltip\" dialog")));
      await session.step(134, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(135, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(136, "And \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(137, "And \"Row Tooltip\" property of grid should be \"target\"", () => propertyShouldBe(page, "Row Tooltip", el("grid"), "target"));
      await session.step(138, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(139, "Then the tooltip should show columns \"target\"", () => tooltipColumns(page, "target"));
      await session.step(140, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(141, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("View > Edit Tooltip leaves a grid set to \"do not show\" as it is (core#474)", async () => {
      await session.step(145, "When user closes all views", () => closeAllViews(page));
      await session.step(146, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
      await session.step(153, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","do not show"]]), [["Show Tooltip","do not show"]]);
      await session.step(155, "And user picks \"View > Edit Tooltip...\" from the top menu", () => pickFromTopMenu(page, "View > Edit Tooltip..."));
      await session.step(156, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(157, "And \"Use tooltip\" input in \"Edit Tooltip\" dialog should not be visible", () => shouldNotBe(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(158, "When user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(159, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(160, "And \"Show Tooltip\" property of grid should be \"do not show\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "do not show"));
      await session.step(161, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(162, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(163, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(164, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(165, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid takes the table's tooltip, and Show Visible Columns In Tooltip is off", async () => {
      await session.step(168, "When user closes all views", () => closeAllViews(page));
      await session.step(169, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
      await session.step(176, "Then \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(177, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(178, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(179, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(180, "Given \"Tooltip\" category in context panel is expanded", () => isExpanded(page, el("\"Tooltip\" category in context panel")));
      await session.step(181, "When user selects \"inherit from table\" in \"Show Tooltip\" property in context panel", () => selectIn(page, "inherit from table", el("\"Show Tooltip\" property in context panel")));
      await session.step(182, "Then \"Show Tooltip\" property of grid should be \"inherit from table\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "inherit from table"));
      await session.step(183, "When user selects \"Always\" in \"Show Column Names\" property in context panel", () => selectIn(page, "Always", el("\"Show Column Names\" property in context panel")));
      await session.step(184, "Then \"Show Column Names\" property of grid should be \"Always\"", () => propertyShouldBe(page, "Show Column Names", el("grid"), "Always"));
      await session.step(185, "And \"Show Visible Columns In Tooltip\" property in context panel should be unchecked", () => shouldBe(page, el("\"Show Visible Columns In Tooltip\" property in context panel"), "unchecked"));
      await session.step(186, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Off, no tooltip while every column fits", async () => {
      await session.step(189, "Then the \"column order\" reading of grid should be \"value, source, target\"", () => readingReads(page, "column order", el("grid"), "value, source, target"));
      await session.step(190, "And grid should have a \"cell 1 of target\" area", () => hasArea(page, el("grid"), "cell 1 of target"));
      await session.step(191, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(192, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(193, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(194, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(195, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched on, the tooltip lists every column while every column fits", async () => {
      await session.step(198, "When user checks \"Show Visible Columns In Tooltip\" property in context panel", () => check(page, el("\"Show Visible Columns In Tooltip\" property in context panel")));
      await session.step(199, "Then \"Show Visible Columns In Tooltip\" property of grid should be \"true\"", () => propertyShouldBe(page, "Show Visible Columns In Tooltip", el("grid"), "true"));
      await session.step(200, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(201, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(202, "And the tooltip should show columns \"value, source, target\"", () => tooltipColumns(page, "value, source, target"));
      await session.step(203, "And the tooltip should show \"source\" as \"Agricultural 'waste'\"", () => tooltipValue(page, "source", "Agricultural 'waste'"));
      await session.step(204, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(205, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched on, the tooltip stays the same with columns pushed off the grid", async () => {
      await session.step(208, "When user remembers the \"column width of value\" reading of grid", () => rememberReading(page, "column width of value", el("grid")));
      await session.step(209, "And user drags the \"column resizer value\" area of grid by 1000 pixels to the right", () => dragAreaBy(page, "column resizer value", el("grid"), 1000, "right"));
      await session.step(210, "Then the \"column width of value\" reading of grid should be higher than remembered", () => readingHigherThanRemembered(page, "column width of value", el("grid")));
      await session.step(211, "And grid should not have a \"cell 1 of source\" area", () => hasNoArea(page, el("grid"), "cell 1 of source"));
      await session.step(212, "And grid should not have a \"cell 1 of target\" area", () => hasNoArea(page, el("grid"), "cell 1 of target"));
      await session.step(213, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(214, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(215, "And the tooltip should show columns \"value, source, target\"", () => tooltipColumns(page, "value, source, target"));
      await session.step(216, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(217, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Off again, the tooltip lists only the columns pushed off the grid", async () => {
      await session.step(220, "When user unchecks \"Show Visible Columns In Tooltip\" property in context panel", () => uncheck(page, el("\"Show Visible Columns In Tooltip\" property in context panel")));
      await session.step(221, "Then \"Show Visible Columns In Tooltip\" property of grid should be \"false\"", () => propertyShouldBe(page, "Show Visible Columns In Tooltip", el("grid"), "false"));
      await session.step(222, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(223, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(224, "And the tooltip should show columns \"source, target\"", () => tooltipColumns(page, "source, target"));
      await session.step(225, "And the tooltip should show \"target\" as \"Bio-conversion\"", () => tooltipValue(page, "target", "Bio-conversion"));
      await session.step(226, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(227, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
