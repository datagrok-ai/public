/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {doubleClickEmptySpace, zoomValueAxis} from '../../../bindings/box-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, closeContextMenu, eventFired, hasArea, hasNoArea, hoverArea, lessInk, listenFor, moreInk, narrowerRange, noErrors, openContextMenu, painted, pickFromContextMenu, pointerAway, propertyShouldBe, rememberRange, rememberedRange, repainted, repaintedBy, resizeTo, resizeWidth, restoreSize, rightClickArea, setProperties, setProperty, tooltipColumns, tooltipNotColumns, tooltipSomeColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot property surface", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot.feature", import.meta.url);
  test("Box plot property surface", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 13, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]));
    await run.scenario("Context menus as property paths", async () => {
      await session.step(19, "Then \"Show Inside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "true"));
      await session.step(20, "When user picks \"Misc > Show Inside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Inside Values", el("box plot viewer")));
      await session.step(21, "Then \"Show Inside Values\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "false"));
      await session.step(22, "And box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await session.step(23, "When user picks \"Misc > Show Outside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Outside Values", el("box plot viewer")));
      await session.step(24, "Then \"Show Outside Values\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Outside Values", el("box plot viewer"), "false"));
      await session.step(25, "And box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await session.step(26, "When user picks \"Misc > Show Inside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Inside Values", el("box plot viewer")));
      await session.step(27, "And user picks \"Misc > Show Outside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Outside Values", el("box plot viewer")));
      await session.step(28, "Then \"Show Inside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "true"));
      await session.step(29, "And \"Show Outside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Outside Values", el("box plot viewer"), "true"));
      await session.step(30, "And box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
      await session.step(31, "When user sets \"Marker Size Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), "WEIGHT"));
      await session.step(32, "And user opens the context menu of box plot viewer", () => openContextMenu(page, el("box plot viewer")));
      await session.step(33, "And user hovers over Markers menu item in context menu", () => hoverOver(page, el("Markers menu item in context menu")));
      await session.step(34, "Then \"Markers > Size\" menu item in context menu should be disabled", () => shouldBe(page, el("\"Markers > Size\" menu item in context menu"), "disabled"));
      await session.step(35, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(36, "And user sets \"Marker Size Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), ""));
    });
    await run.scenario("Statistics and group-comparison menu regions", async () => {
      await session.step(39, "Given user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Statistics","true"],["Show Group Comparison","false"],["Show P Value","true"]]));
      await session.step(43, "When user right-clicks on the \"stats\" area of box plot viewer", () => rightClickArea(page, "stats", el("box plot viewer")));
      await session.step(44, "And user hovers over \"Group Comparison\" menu item in context menu", () => hoverOver(page, el("\"Group Comparison\" menu item in context menu")));
      await session.step(45, "Then \"Show Assumption Checks\" menu item in context menu should be disabled", () => shouldBe(page, el("\"Show Assumption Checks\" menu item in context menu"), "disabled"));
      await session.step(46, "When user hovers over \"Show Assumption Checks\" menu item in context menu", () => hoverOver(page, el("\"Show Assumption Checks\" menu item in context menu")));
      await session.step(47, "Then tooltip should contain text \"Show Group Comparison\"", () => shouldContainText(page, el("tooltip"), "Show Group Comparison"));
      await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(49, "And user right-clicks on the \"p value\" area of box plot viewer", () => rightClickArea(page, "p value", el("box plot viewer")));
      await session.step(50, "Then context menu should contain text \"Show P Value\"", () => shouldContainText(page, el("context menu"), "Show P Value"));
      await session.step(51, "And context menu should not contain text \"Statistics Format\"", () => shouldNotContainText(page, el("context menu"), "Statistics Format"));
      await session.step(52, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(53, "And user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(54, "And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(55, "And user clicks on show-group-stats icon in box plot viewer", () => clickOn(page, el("show-group-stats icon in box plot viewer")));
      await session.step(56, "Then \"Show Group Comparison\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Group Comparison", el("box plot viewer"), "true"));
      await session.step(57, "When user right-clicks on the \"group comparison\" area of box plot viewer", () => rightClickArea(page, "group comparison", el("box plot viewer")));
      await session.step(58, "Then \"Add t-Test Table\" menu item in context menu should be visible", () => shouldBe(page, el("\"Add t-Test Table\" menu item in context menu"), "visible"));
      await session.step(59, "And \"Show Assumption Checks\" menu item in context menu should be enabled", () => shouldBe(page, el("\"Show Assumption Checks\" menu item in context menu"), "enabled"));
      await session.step(60, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(61, "And user sets \"Show Group Comparison\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Group Comparison", el("box plot viewer"), "false"));
    });
    await run.scenario("Resize and auto layout", async () => {
      await session.step(64, "Given user sets \"Auto Layout\" property of box plot viewer to \"true\"", () => setProperty(page, "Auto Layout", el("box plot viewer"), "true"));
      await session.step(65, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(66, "Then \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
      await session.step(67, "When user resizes box plot viewer to 170 by 150", () => resizeTo(page, el("box plot viewer"), 170, 150));
      await session.step(68, "Then \"Marker Color\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "hidden"));
      await session.step(69, "When user restores the size of box plot viewer", () => restoreSize(page, el("box plot viewer")));
      await session.step(70, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(71, "Then \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
      await session.step(72, "When user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(73, "And user resizes box plot viewer to 120 wide", () => resizeWidth(page, el("box plot viewer"), 120));
      await session.step(74, "And user restores the size of box plot viewer", () => restoreSize(page, el("box plot viewer")));
      await session.step(75, "Then no errors should have been logged", () => noErrors(page));
      await session.step(76, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Marker gate and size scaling", async () => {
      await session.step(79, "When user sets \"Show Markers\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Markers", el("box plot viewer"), "false"));
      await session.step(80, "Then box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await session.step(81, "When user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
      await session.step(82, "Then \"Marker Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "disabled"));
      await session.step(83, "When user sets \"Show Markers\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Markers", el("box plot viewer"), "true"));
      await session.step(84, "Then box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
      await session.step(85, "And \"Marker Type\" property in context panel should be enabled", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "enabled"));
      await session.step(86, "When user sets \"Marker Size Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), "WEIGHT"));
      await session.step(87, "And user sets \"Marker Size Scaling\" property of box plot viewer to \"logarithmic\"", () => setProperty(page, "Marker Size Scaling", el("box plot viewer"), "logarithmic"));
      await session.step(88, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(89, "When user sets \"Marker Size Scaling\" property of box plot viewer to \"linear\"", () => setProperty(page, "Marker Size Scaling", el("box plot viewer"), "linear"));
      await session.step(90, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(91, "When user sets \"Marker Size Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), ""));
    });
    await run.scenario("Whisker and control-band style", async () => {
      await session.step(94, "When user sets \"Whisker Line Width\" property of box plot viewer to \"4\"", () => setProperty(page, "Whisker Line Width", el("box plot viewer"), "4"));
      await session.step(95, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(96, "When user sets \"Whisker Width Ratio\" property of box plot viewer to \"0.3\"", () => setProperty(page, "Whisker Width Ratio", el("box plot viewer"), "0.3"));
      await session.step(97, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(98, "When user sets \"Control Band Color\" property of box plot viewer to \"#00AA00\"", () => setProperty(page, "Control Band Color", el("box plot viewer"), "#00AA00"));
      await session.step(99, "Then no errors should have been logged", () => noErrors(page));
      await session.step(100, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Whisker Line Width","2"],["Whisker Width Ratio","0.5"]]));
    });
    await run.scenario("Controls visibility", async () => {
      await session.step(105, "Then \"Show Size Selector\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Size Selector", el("box plot viewer"), "false"));
      await session.step(106, "And \"Marker Size\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "hidden"));
      await session.step(107, "When user sets \"Show Size Selector\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Size Selector", el("box plot viewer"), "true"));
      await session.step(108, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(109, "Then \"Marker Size\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "visible"));
      await session.step(110, "When user sets \"Show Value Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Value Selector", el("box plot viewer"), "false"));
      await session.step(111, "Then Value column input in box plot viewer should be hidden", () => shouldBe(page, el("Value column input in box plot viewer"), "hidden"));
      await session.step(112, "When user sets \"Show Color Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Color Selector", el("box plot viewer"), "false"));
      await session.step(113, "Then \"Marker Color\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "hidden"));
      await session.step(114, "When user sets \"Show Category Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Category Selector", el("box plot viewer"), "false"));
      await session.step(115, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(116, "When user sets \"Show Value Axis\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Value Axis", el("box plot viewer"), "false"));
      await session.step(117, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(118, "And box plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("box plot viewer"), "y axis"));
      await session.step(119, "When user sets \"Show Category Axis\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Category Axis", el("box plot viewer"), "false"));
      await session.step(120, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(121, "And box plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("box plot viewer"), "x axis"));
      await session.step(122, "When user sets \"Show Category Axis\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Category Axis", el("box plot viewer"), "true"));
      await session.step(123, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(124, "And box plot viewer should have an \"x axis\" area", () => hasArea(page, el("box plot viewer"), "x axis"));
      await session.step(125, "When user sets \"Show Value Axis\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Value Axis", el("box plot viewer"), "true"));
      await session.step(126, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(127, "And box plot viewer should have a \"y axis\" area", () => hasArea(page, el("box plot viewer"), "y axis"));
      await session.step(128, "When user sets \"Show Category Selector\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Category Selector", el("box plot viewer"), "true"));
      await session.step(129, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(130, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Value Selector","true"],["Show Color Selector","true"],["Show Size Selector","false"]]));
      await session.step(134, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(135, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(136, "Then \"Marker Size\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "hidden"));
      await session.step(137, "And Value column input in box plot viewer should be visible", () => shouldBe(page, el("Value column input in box plot viewer"), "visible"));
      await session.step(138, "And \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
    });
    await run.scenario("Title and description", async () => {
      await session.step(141, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Title","true"],["Title","Age by Race"]]));
      await session.step(144, "Then title of box plot viewer should have text \"Age by Race\"", () => shouldHaveText(page, el("title of box plot viewer"), "Age by Race"));
      await session.step(145, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Description","Box plot of patient ages"],["Description Visibility Mode","Always"]]));
      await session.step(148, "Then description of box plot viewer should have text \"Box plot of patient ages\"", () => shouldHaveText(page, el("description of box plot viewer"), "Box plot of patient ages"));
      await session.step(149, "When user sets \"Description Position\" property of box plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("box plot viewer"), "Bottom"));
      await session.step(150, "Then description of box plot viewer should be visible", () => shouldBe(page, el("description of box plot viewer"), "visible"));
      await session.step(151, "When user sets \"Description Visibility Mode\" property of box plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("box plot viewer"), "Never"));
      await session.step(152, "Then description of box plot viewer should be absent", () => shouldBe(page, el("description of box plot viewer"), "absent"));
      await session.step(153, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
    });
    await run.scenario("Axis font", async () => {
      await session.step(161, "When user sets \"Axis Font\" property of box plot viewer to \"normal normal 16px \\\"Roboto\\\"\"", () => setProperty(page, "Axis Font", el("box plot viewer"), "normal normal 16px \"Roboto\""));
      await session.step(162, "Then \"Axis Font\" property of box plot viewer should be \"normal normal 16px \\\"Roboto\\\"\"", () => propertyShouldBe(page, "Axis Font", el("box plot viewer"), "normal normal 16px \"Roboto\""));
      await session.step(163, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(164, "When user sets \"Axis Font\" property of box plot viewer to \"normal normal 10px \\\"Roboto\\\"\"", () => setProperty(page, "Axis Font", el("box plot viewer"), "normal normal 10px \"Roboto\""));
      await session.step(165, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Date category mapping", async () => {
      await session.step(168, "When user sets \"Category 1\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Category 1", el("box plot viewer"), "STARTED"));
      await session.step(169, "Then \"Category 1\" property of box plot viewer should be \"STARTED\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "STARTED"));
      await session.step(170, "When user sets \"Category 1 Map\" property of box plot viewer to \"month\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "month"));
      await session.step(171, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(172, "And \"Category 1 Map\" property of box plot viewer should be \"month\"", () => propertyShouldBe(page, "Category 1 Map", el("box plot viewer"), "month"));
      await session.step(173, "When user sets \"Category 1 Map\" property of box plot viewer to \"quarter\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "quarter"));
      await session.step(174, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(175, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(176, "Then \"Category 1\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(177, "And no errors should have been logged", () => noErrors(page));
      await session.step(178, "When user sets \"Category 1\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Category 1", el("box plot viewer"), "SEX"));
    });
    await run.scenario("Custom tooltip", async () => {
      await session.step(181, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","RACE"],["Marker Size","10"],["Row Tooltip","AGE\nSEX\nWEIGHT"],["Show Tooltip","show custom tooltip"]]));
      await session.step(186, "Then \"Row Tooltip\" property of box plot viewer should be \"AGE\\nSEX\\nWEIGHT\"", () => propertyShouldBe(page, "Row Tooltip", el("box plot viewer"), "AGE\\nSEX\\nWEIGHT"));
      await session.step(187, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(188, "Then the tooltip should show columns \"AGE, SEX, WEIGHT\"", () => tooltipColumns(page, "AGE, SEX, WEIGHT"));
      await session.step(189, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(190, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(191, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Category 1","SEX"]]));
      await session.step(195, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(196, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(197, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(198, "And the tooltip should not show columns \"AGE, SEX, WEIGHT\"", () => tooltipNotColumns(page, "AGE, SEX, WEIGHT"));
      await session.step(199, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
    });
    await run.scenario("Table switching resets Category 2", async () => {
      await session.step(202, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(203, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(204, "When user sets \"Category 2\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 2", el("box plot viewer"), "RACE"));
      await session.step(205, "Then \"Category 2\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Category 2", el("box plot viewer"), "RACE"));
      await session.step(206, "When user sets \"Table\" property of box plot viewer to \"spgi-100\"", () => setProperty(page, "Table", el("box plot viewer"), "spgi-100"));
      await session.step(207, "Then \"Category 2\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Category 2", el("box plot viewer"), ""));
      await session.step(208, "And no errors should have been logged", () => noErrors(page));
      await session.step(209, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","Average Mass"],["Category 1","Series"]]));
      await session.step(212, "Then \"Table\" property of box plot viewer should be \"spgi-100\"", () => propertyShouldBe(page, "Table", el("box plot viewer"), "spgi-100"));
      await session.step(213, "And box plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("box plot viewer"), "spgi-100"));
      await session.step(214, "And box plot viewer should have a \"category Triazoles\" area", () => hasArea(page, el("box plot viewer"), "category Triazoles"));
      await session.step(215, "When user sets \"Table\" property of box plot viewer to \"demog-1000\"", () => setProperty(page, "Table", el("box plot viewer"), "demog-1000"));
      await session.step(216, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","AGE"],["Category 1","SEX"],["Category 2",""]]));
      await session.step(220, "Then \"Table\" property of box plot viewer should be \"demog-1000\"", () => propertyShouldBe(page, "Table", el("box plot viewer"), "demog-1000"));
      await session.step(221, "And box plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("box plot viewer"), "demog-1000"));
    });
    await run.scenario("Coloring keeps the render valid", async () => {
      await session.step(224, "When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(225, "Then \"Marker Color Column\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(226, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(227, "And no errors should have been logged", () => noErrors(page));
      await session.step(228, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Double-click resets the view", async () => {
      await session.step(231, "Given user listens for \"d4-boxplot-reset-view\" event on box plot viewer", () => listenFor(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(232, "And user remembers the value range of box plot viewer", () => rememberRange(page, el("box plot viewer")));
      await session.step(233, "When user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
      await session.step(234, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(235, "When user double-clicks on empty plot space of box plot viewer", () => doubleClickEmptySpace(page));
      await session.step(236, "Then \"d4-boxplot-reset-view\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(237, "And box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
    });
    run.finish();
  });
});
