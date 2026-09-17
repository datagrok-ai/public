/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, closeContextMenu, eventFired, hasArea, hasNoArea, hoverArea, lessInk, listenFor, moreInk, narrowerRange, noErrors, openContextMenu, painted, pickFromContextMenu, pointerAway, propertyShouldBe, rememberRange, rememberedRange, repainted, repaintedBy, resizeTo, resizeWidth, restoreSize, rightClickArea, setProperties, setProperty, tooltipColumns, tooltipNotColumns, tooltipSomeColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doubleClickEmptySpace, zoomValueAxis} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot property surface", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot.feature", import.meta.url);
  test("Box plot property surface", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]));
    await run.scenario("Context menus as property paths", async () => {
      await session.step(20, "Then \"Show Inside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "true"));
      await session.step(21, "When user picks \"Misc > Show Inside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Inside Values", el("box plot viewer")));
      await session.step(22, "Then \"Show Inside Values\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "false"));
      await session.step(23, "And box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await session.step(24, "When user picks \"Misc > Show Outside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Outside Values", el("box plot viewer")));
      await session.step(25, "Then \"Show Outside Values\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Outside Values", el("box plot viewer"), "false"));
      await session.step(26, "And box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await session.step(27, "When user picks \"Misc > Show Inside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Inside Values", el("box plot viewer")));
      await session.step(28, "And user picks \"Misc > Show Outside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Outside Values", el("box plot viewer")));
      await session.step(29, "Then \"Show Inside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "true"));
      await session.step(30, "And \"Show Outside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Outside Values", el("box plot viewer"), "true"));
      await session.step(31, "And box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
      await session.step(32, "When user sets \"Marker Size Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), "WEIGHT"));
      await session.step(33, "And user opens the context menu of box plot viewer", () => openContextMenu(page, el("box plot viewer")));
      await session.step(34, "And user hovers over Markers menu item in context menu", () => hoverOver(page, el("Markers menu item in context menu")));
      await session.step(35, "Then \"Markers > Size\" menu item in context menu should be disabled", () => shouldBe(page, el("\"Markers > Size\" menu item in context menu"), "disabled"));
      await session.step(36, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(37, "And user sets \"Marker Size Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), ""));
    });
    await run.scenario("Statistics and group-comparison menu regions", async () => {
      await session.step(40, "Given user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Statistics","true"],["Show Group Comparison","false"],["Show P Value","true"]]));
      await session.step(44, "When user right-clicks on the \"stats\" area of box plot viewer", () => rightClickArea(page, "stats", el("box plot viewer")));
      await session.step(45, "And user hovers over \"Group Comparison\" menu item in context menu", () => hoverOver(page, el("\"Group Comparison\" menu item in context menu")));
      await session.step(46, "Then \"Show Assumption Checks\" menu item in context menu should be disabled", () => shouldBe(page, el("\"Show Assumption Checks\" menu item in context menu"), "disabled"));
      await session.step(47, "When user hovers over \"Show Assumption Checks\" menu item in context menu", () => hoverOver(page, el("\"Show Assumption Checks\" menu item in context menu")));
      await session.step(48, "Then tooltip should contain text \"Show Group Comparison\"", () => shouldContainText(page, el("tooltip"), "Show Group Comparison"));
      await session.step(49, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(50, "And user right-clicks on the \"p value\" area of box plot viewer", () => rightClickArea(page, "p value", el("box plot viewer")));
      await session.step(51, "Then context menu should contain text \"Show P Value\"", () => shouldContainText(page, el("context menu"), "Show P Value"));
      await session.step(52, "And context menu should not contain text \"Statistics Format\"", () => shouldNotContainText(page, el("context menu"), "Statistics Format"));
      await session.step(53, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(54, "And user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(55, "And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(56, "And user clicks on show-group-stats icon in box plot viewer", () => clickOn(page, el("show-group-stats icon in box plot viewer")));
      await session.step(57, "Then \"Show Group Comparison\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Group Comparison", el("box plot viewer"), "true"));
      await session.step(58, "When user right-clicks on the \"group comparison\" area of box plot viewer", () => rightClickArea(page, "group comparison", el("box plot viewer")));
      await session.step(59, "Then \"Add t-Test Table\" menu item in context menu should be visible", () => shouldBe(page, el("\"Add t-Test Table\" menu item in context menu"), "visible"));
      await session.step(60, "And \"Show Assumption Checks\" menu item in context menu should be enabled", () => shouldBe(page, el("\"Show Assumption Checks\" menu item in context menu"), "enabled"));
      await session.step(61, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(62, "And user sets \"Show Group Comparison\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Group Comparison", el("box plot viewer"), "false"));
    });
    await run.scenario("Resize and auto layout", async () => {
      await session.step(65, "Given user sets \"Auto Layout\" property of box plot viewer to \"true\"", () => setProperty(page, "Auto Layout", el("box plot viewer"), "true"));
      await session.step(66, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(67, "Then \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
      await session.step(68, "When user resizes box plot viewer to 170 by 150", () => resizeTo(page, el("box plot viewer"), 170, 150));
      await session.step(69, "Then \"Marker Color\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "hidden"));
      await session.step(70, "When user restores the size of box plot viewer", () => restoreSize(page, el("box plot viewer")));
      await session.step(71, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(72, "Then \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
      await session.step(73, "When user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(74, "And user resizes box plot viewer to 120 wide", () => resizeWidth(page, el("box plot viewer"), 120));
      await session.step(75, "And user restores the size of box plot viewer", () => restoreSize(page, el("box plot viewer")));
      await session.step(76, "Then no errors should have been logged", () => noErrors(page));
      await session.step(77, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Marker gate and size scaling", async () => {
      await session.step(80, "When user sets \"Show Markers\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Markers", el("box plot viewer"), "false"));
      await session.step(81, "Then box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await session.step(82, "When user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
      await session.step(83, "Then \"Marker Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "disabled"));
      await session.step(84, "When user sets \"Show Markers\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Markers", el("box plot viewer"), "true"));
      await session.step(85, "Then box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
      await session.step(86, "And \"Marker Type\" property in context panel should be enabled", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "enabled"));
      await session.step(87, "When user sets \"Marker Size Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), "WEIGHT"));
      await session.step(88, "And user sets \"Marker Size Scaling\" property of box plot viewer to \"logarithmic\"", () => setProperty(page, "Marker Size Scaling", el("box plot viewer"), "logarithmic"));
      await session.step(89, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(90, "When user sets \"Marker Size Scaling\" property of box plot viewer to \"linear\"", () => setProperty(page, "Marker Size Scaling", el("box plot viewer"), "linear"));
      await session.step(91, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(92, "When user sets \"Marker Size Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), ""));
    });
    await run.scenario("Whisker and control-band style", async () => {
      await session.step(95, "When user sets \"Whisker Line Width\" property of box plot viewer to \"4\"", () => setProperty(page, "Whisker Line Width", el("box plot viewer"), "4"));
      await session.step(96, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(97, "When user sets \"Whisker Width Ratio\" property of box plot viewer to \"0.3\"", () => setProperty(page, "Whisker Width Ratio", el("box plot viewer"), "0.3"));
      await session.step(98, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(99, "When user sets \"Control Band Color\" property of box plot viewer to \"#00AA00\"", () => setProperty(page, "Control Band Color", el("box plot viewer"), "#00AA00"));
      await session.step(100, "Then no errors should have been logged", () => noErrors(page));
      await session.step(101, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Whisker Line Width","2"],["Whisker Width Ratio","0.5"]]));
    });
    await run.scenario("Controls visibility", async () => {
      await session.step(106, "Then \"Show Size Selector\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Size Selector", el("box plot viewer"), "false"));
      await session.step(107, "And \"Marker Size\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "hidden"));
      await session.step(108, "When user sets \"Show Size Selector\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Size Selector", el("box plot viewer"), "true"));
      await session.step(109, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(110, "Then \"Marker Size\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "visible"));
      await session.step(111, "When user sets \"Show Value Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Value Selector", el("box plot viewer"), "false"));
      await session.step(112, "Then Value column input in box plot viewer should be hidden", () => shouldBe(page, el("Value column input in box plot viewer"), "hidden"));
      await session.step(113, "When user sets \"Show Color Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Color Selector", el("box plot viewer"), "false"));
      await session.step(114, "Then \"Marker Color\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "hidden"));
      await session.step(115, "When user sets \"Show Category Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Category Selector", el("box plot viewer"), "false"));
      await session.step(116, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(117, "When user sets \"Show Value Axis\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Value Axis", el("box plot viewer"), "false"));
      await session.step(118, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(119, "And box plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("box plot viewer"), "y axis"));
      await session.step(120, "When user sets \"Show Category Axis\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Category Axis", el("box plot viewer"), "false"));
      await session.step(121, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(122, "And box plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("box plot viewer"), "x axis"));
      await session.step(123, "When user sets \"Show Category Axis\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Category Axis", el("box plot viewer"), "true"));
      await session.step(124, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(125, "And box plot viewer should have an \"x axis\" area", () => hasArea(page, el("box plot viewer"), "x axis"));
      await session.step(126, "When user sets \"Show Value Axis\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Value Axis", el("box plot viewer"), "true"));
      await session.step(127, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(128, "And box plot viewer should have a \"y axis\" area", () => hasArea(page, el("box plot viewer"), "y axis"));
      await session.step(129, "When user sets \"Show Category Selector\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Category Selector", el("box plot viewer"), "true"));
      await session.step(130, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(131, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Value Selector","true"],["Show Color Selector","true"],["Show Size Selector","false"]]));
      await session.step(135, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(136, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(137, "Then \"Marker Size\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "hidden"));
      await session.step(138, "And Value column input in box plot viewer should be visible", () => shouldBe(page, el("Value column input in box plot viewer"), "visible"));
      await session.step(139, "And \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
    });
    await run.scenario("Axis font", async () => {
      await session.step(142, "When user sets \"Axis Font\" property of box plot viewer to \"normal normal 16px \\\"Roboto\\\"\"", () => setProperty(page, "Axis Font", el("box plot viewer"), "normal normal 16px \"Roboto\""));
      await session.step(143, "Then \"Axis Font\" property of box plot viewer should be \"normal normal 16px \\\"Roboto\\\"\"", () => propertyShouldBe(page, "Axis Font", el("box plot viewer"), "normal normal 16px \"Roboto\""));
      await session.step(144, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(145, "When user sets \"Axis Font\" property of box plot viewer to \"normal normal 10px \\\"Roboto\\\"\"", () => setProperty(page, "Axis Font", el("box plot viewer"), "normal normal 10px \"Roboto\""));
      await session.step(146, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Date category mapping", async () => {
      await session.step(149, "When user sets \"Category 1\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Category 1", el("box plot viewer"), "STARTED"));
      await session.step(150, "Then \"Category 1\" property of box plot viewer should be \"STARTED\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "STARTED"));
      await session.step(151, "When user sets \"Category 1 Map\" property of box plot viewer to \"month\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "month"));
      await session.step(152, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(153, "And \"Category 1 Map\" property of box plot viewer should be \"month\"", () => propertyShouldBe(page, "Category 1 Map", el("box plot viewer"), "month"));
      await session.step(154, "When user sets \"Category 1 Map\" property of box plot viewer to \"quarter\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "quarter"));
      await session.step(155, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(156, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(157, "Then \"Category 1\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
      await session.step(159, "When user sets \"Category 1\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Category 1", el("box plot viewer"), "SEX"));
    });
    await run.scenario("Custom tooltip", async () => {
      await session.step(162, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","RACE"],["Marker Size","10"],["Row Tooltip","AGE\nSEX\nWEIGHT"],["Show Tooltip","show custom tooltip"]]));
      await session.step(167, "Then \"Row Tooltip\" property of box plot viewer should be \"AGE\\nSEX\\nWEIGHT\"", () => propertyShouldBe(page, "Row Tooltip", el("box plot viewer"), "AGE\\nSEX\\nWEIGHT"));
      await session.step(168, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(169, "Then the tooltip should show columns \"AGE, SEX, WEIGHT\"", () => tooltipColumns(page, "AGE, SEX, WEIGHT"));
      await session.step(170, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(171, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(172, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Category 1","SEX"]]));
      await session.step(176, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(177, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(178, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(179, "And the tooltip should not show columns \"AGE, SEX, WEIGHT\"", () => tooltipNotColumns(page, "AGE, SEX, WEIGHT"));
      await session.step(180, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
    });
    await run.scenario("Table switching resets Category 2", async () => {
      await session.step(183, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(184, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(185, "When user sets \"Category 2\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 2", el("box plot viewer"), "RACE"));
      await session.step(186, "Then \"Category 2\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Category 2", el("box plot viewer"), "RACE"));
      await session.step(187, "When user sets \"Table\" property of box plot viewer to \"spgi-100\"", () => setProperty(page, "Table", el("box plot viewer"), "spgi-100"));
      await session.step(188, "Then \"Category 2\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Category 2", el("box plot viewer"), ""));
      await session.step(189, "And no errors should have been logged", () => noErrors(page));
      await session.step(190, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","Average Mass"],["Category 1","Series"]]));
      await session.step(193, "Then \"Table\" property of box plot viewer should be \"spgi-100\"", () => propertyShouldBe(page, "Table", el("box plot viewer"), "spgi-100"));
      await session.step(194, "And box plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("box plot viewer"), "spgi-100"));
      await session.step(195, "And box plot viewer should have a \"category Triazoles\" area", () => hasArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(196, "When user sets \"Table\" property of box plot viewer to \"demog-1000\"", () => setProperty(page, "Table", el("box plot viewer"), "demog-1000"));
      await session.step(197, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","AGE"],["Category 1","SEX"],["Category 2",""]]));
      await session.step(201, "Then \"Table\" property of box plot viewer should be \"demog-1000\"", () => propertyShouldBe(page, "Table", el("box plot viewer"), "demog-1000"));
      await session.step(202, "And box plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("box plot viewer"), "demog-1000"));
    });
    await run.scenario("Coloring keeps the render valid", async () => {
      await session.step(205, "When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(206, "Then \"Marker Color Column\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(207, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(208, "And no errors should have been logged", () => noErrors(page));
      await session.step(209, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Double-click resets the view", async () => {
      await session.step(212, "Given user listens for \"d4-boxplot-reset-view\" event on box plot viewer", () => listenFor(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(213, "And user remembers the value range of box plot viewer", () => rememberRange(page, el("box plot viewer")));
      await session.step(214, "When user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
      await session.step(215, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(216, "When user double-clicks on empty plot space of box plot viewer", () => doubleClickEmptySpace(page, el("box plot viewer")));
      await session.step(217, "Then \"d4-boxplot-reset-view\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(218, "And box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
    });
    run.finish();
  });
});
