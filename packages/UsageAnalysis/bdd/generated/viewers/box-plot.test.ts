/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {doubleClickEmptySpace, fullRange, narrowedRange, zoomValueAxis} from '../../bindings/box-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, eventFired, hoverArea, lessInk, listenFor, moreInk, noErrors, openContextMenu, painted, pickFromContextMenu, pointerAway, propertyShouldBe, propertyShouldNotBe, repainted, resizeTo, resizeWidth, restoreSize, rightClickArea, setProperties, setProperty, tooltipColumns, tooltipNotColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot property surface", () => {
  const session = feature(test);
  test("Box plot property surface", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 13);
    await test.step("Given user is logged in", () => loggedIn(page));
    await test.step("And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await test.step("And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]));
    await run.scenario("Context menus as property paths", async () => {
      await test.step("Then \"Show Inside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "true"));
      await test.step("When user picks \"Misc > Show Inside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Inside Values", el("box plot viewer")));
      await test.step("Then \"Show Inside Values\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "false"));
      await test.step("And box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await test.step("When user picks \"Misc > Show Outside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Outside Values", el("box plot viewer")));
      await test.step("Then \"Show Outside Values\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Outside Values", el("box plot viewer"), "false"));
      await test.step("And box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await test.step("When user picks \"Misc > Show Inside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Inside Values", el("box plot viewer")));
      await test.step("And user picks \"Misc > Show Outside Values\" from the context menu of box plot viewer", () => pickFromContextMenu(page, "Misc > Show Outside Values", el("box plot viewer")));
      await test.step("Then \"Show Inside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Inside Values", el("box plot viewer"), "true"));
      await test.step("And \"Show Outside Values\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Outside Values", el("box plot viewer"), "true"));
      await test.step("And box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
      await test.step("When user sets \"Marker Size Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), "WEIGHT"));
      await test.step("And user opens the context menu of box plot viewer", () => openContextMenu(page, el("box plot viewer")));
      await test.step("And user hovers over Markers menu item in context menu", () => hoverOver(page, el("Markers menu item in context menu")));
      await test.step("Then \"Markers > Size\" menu item in context menu should be disabled", () => shouldBe(page, el("\"Markers > Size\" menu item in context menu"), "disabled"));
      await test.step("When user closes the context menu", () => closeContextMenu(page));
      await test.step("And user sets \"Marker Size Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), ""));
    });
    await run.scenario("Statistics and group-comparison menu regions", async () => {
      await test.step("Given user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Statistics","true"],["Show Group Comparison","false"],["Show P Value","true"]]));
      await test.step("When user right-clicks on the \"stats\" area of box plot viewer", () => rightClickArea(page, "stats", el("box plot viewer")));
      await test.step("And user hovers over \"Group Comparison\" menu item in context menu", () => hoverOver(page, el("\"Group Comparison\" menu item in context menu")));
      await test.step("Then \"Show Assumption Checks\" menu item in context menu should be disabled", () => shouldBe(page, el("\"Show Assumption Checks\" menu item in context menu"), "disabled"));
      await test.step("When user hovers over \"Show Assumption Checks\" menu item in context menu", () => hoverOver(page, el("\"Show Assumption Checks\" menu item in context menu")));
      await test.step("Then tooltip should contain text \"Show Group Comparison\"", () => shouldContainText(page, el("tooltip"), "Show Group Comparison"));
      await test.step("When user closes the context menu", () => closeContextMenu(page));
      await test.step("And user right-clicks on the \"p value\" area of box plot viewer", () => rightClickArea(page, "p value", el("box plot viewer")));
      await test.step("Then context menu should contain text \"Show P Value\"", () => shouldContainText(page, el("context menu"), "Show P Value"));
      await test.step("And context menu should not contain text \"Statistics Format\"", () => shouldNotContainText(page, el("context menu"), "Statistics Format"));
      await test.step("When user closes the context menu", () => closeContextMenu(page));
      await test.step("And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await test.step("And user clicks on show-group-stats icon in box plot viewer", () => clickOn(page, el("show-group-stats icon in box plot viewer")));
      await test.step("Then \"Show Group Comparison\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Group Comparison", el("box plot viewer"), "true"));
      await test.step("When user right-clicks on the \"group comparison\" area of box plot viewer", () => rightClickArea(page, "group comparison", el("box plot viewer")));
      await test.step("Then context menu should contain text \"Table\"", () => shouldContainText(page, el("context menu"), "Table"));
      await test.step("And \"Show Assumption Checks\" menu item in context menu should be enabled", () => shouldBe(page, el("\"Show Assumption Checks\" menu item in context menu"), "enabled"));
      await test.step("When user closes the context menu", () => closeContextMenu(page));
      await test.step("And user sets \"Show Group Comparison\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Group Comparison", el("box plot viewer"), "false"));
    });
    await run.scenario("Resize and auto layout", async () => {
      await test.step("Given user sets \"Auto Layout\" property of box plot viewer to \"true\"", () => setProperty(page, "Auto Layout", el("box plot viewer"), "true"));
      await test.step("When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await test.step("Then \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
      await test.step("When user resizes box plot viewer to 170 by 150", () => resizeTo(page, el("box plot viewer"), 170, 150));
      await test.step("Then \"Marker Color\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "hidden"));
      await test.step("When user restores the size of box plot viewer", () => restoreSize(page, el("box plot viewer")));
      await test.step("And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await test.step("Then \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
      await test.step("When user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await test.step("And user resizes box plot viewer to 120 wide", () => resizeWidth(page, el("box plot viewer"), 120));
      await test.step("And user restores the size of box plot viewer", () => restoreSize(page, el("box plot viewer")));
      await test.step("Then no errors should have been logged", () => noErrors(page));
      await test.step("When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Marker gate and size scaling", async () => {
      await test.step("When user sets \"Show Markers\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Markers", el("box plot viewer"), "false"));
      await test.step("Then box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
      await test.step("When user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
      await test.step("Then \"Marker Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "disabled"));
      await test.step("When user sets \"Show Markers\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Markers", el("box plot viewer"), "true"));
      await test.step("Then box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
      await test.step("And \"Marker Type\" property in context panel should be enabled", () => shouldBe(page, el("\"Marker Type\" property in context panel"), "enabled"));
      await test.step("When user sets \"Marker Size Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), "WEIGHT"));
      await test.step("And user sets \"Marker Size Scaling\" property of box plot viewer to \"logarithmic\"", () => setProperty(page, "Marker Size Scaling", el("box plot viewer"), "logarithmic"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Marker Size Scaling\" property of box plot viewer to \"linear\"", () => setProperty(page, "Marker Size Scaling", el("box plot viewer"), "linear"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Marker Size Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Size Column", el("box plot viewer"), ""));
    });
    await run.scenario("Whisker and control-band style", async () => {
      await test.step("When user sets \"Whisker Line Width\" property of box plot viewer to \"4\"", () => setProperty(page, "Whisker Line Width", el("box plot viewer"), "4"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Whisker Width Ratio\" property of box plot viewer to \"0.3\"", () => setProperty(page, "Whisker Width Ratio", el("box plot viewer"), "0.3"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Control Band Color\" property of box plot viewer to \"#00AA00\"", () => setProperty(page, "Control Band Color", el("box plot viewer"), "#00AA00"));
      await test.step("Then no errors should have been logged", () => noErrors(page));
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Whisker Line Width","2"],["Whisker Width Ratio","0.5"]]));
    });
    await run.scenario("Controls visibility", async () => {
      await test.step("Then \"Show Size Selector\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Size Selector", el("box plot viewer"), "false"));
      await test.step("And \"Marker Size\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "hidden"));
      await test.step("When user sets \"Show Size Selector\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Size Selector", el("box plot viewer"), "true"));
      await test.step("And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await test.step("Then \"Marker Size\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "visible"));
      await test.step("When user sets \"Show Value Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Value Selector", el("box plot viewer"), "false"));
      await test.step("Then Value column input in box plot viewer should be hidden", () => shouldBe(page, el("Value column input in box plot viewer"), "hidden"));
      await test.step("When user sets \"Show Color Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Color Selector", el("box plot viewer"), "false"));
      await test.step("Then \"Marker Color\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "hidden"));
      await test.step("When user sets \"Show Category Selector\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Category Selector", el("box plot viewer"), "false"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Show Value Axis\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Value Axis", el("box plot viewer"), "false"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Show Category Axis\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Category Axis", el("box plot viewer"), "false"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Category Axis","true"],["Show Value Axis","true"],["Show Category Selector","true"],["Show Value Selector","true"],["Show Color Selector","true"],["Show Size Selector","false"]]));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await test.step("Then \"Marker Size\" column input in box plot viewer should be hidden", () => shouldBe(page, el("\"Marker Size\" column input in box plot viewer"), "hidden"));
      await test.step("And Value column input in box plot viewer should be visible", () => shouldBe(page, el("Value column input in box plot viewer"), "visible"));
      await test.step("And \"Marker Color\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Marker Color\" column input in box plot viewer"), "visible"));
    });
    await run.scenario("Title and description", async () => {
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Title","true"],["Title","Age by Race"]]));
      await test.step("Then title of box plot viewer should have text \"Age by Race\"", () => shouldHaveText(page, el("title of box plot viewer"), "Age by Race"));
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Description","Box plot of patient ages"],["Description Visibility Mode","Always"]]));
      await test.step("Then description of box plot viewer should have text \"Box plot of patient ages\"", () => shouldHaveText(page, el("description of box plot viewer"), "Box plot of patient ages"));
      await test.step("When user sets \"Description Position\" property of box plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("box plot viewer"), "Bottom"));
      await test.step("Then description of box plot viewer should be visible", () => shouldBe(page, el("description of box plot viewer"), "visible"));
      await test.step("When user sets \"Description Visibility Mode\" property of box plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("box plot viewer"), "Never"));
      await test.step("Then description of box plot viewer should be absent", () => shouldBe(page, el("description of box plot viewer"), "absent"));
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
    });
    await run.scenario("Axis font", async () => {
      await test.step("When user sets \"Axis Font\" property of box plot viewer to \"normal normal 16px \\\"Roboto\\\"\"", () => setProperty(page, "Axis Font", el("box plot viewer"), "normal normal 16px \"Roboto\""));
      await test.step("Then \"Axis Font\" property of box plot viewer should be \"normal normal 16px \\\"Roboto\\\"\"", () => propertyShouldBe(page, "Axis Font", el("box plot viewer"), "normal normal 16px \"Roboto\""));
      await test.step("And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Axis Font\" property of box plot viewer to \"normal normal 10px \\\"Roboto\\\"\"", () => setProperty(page, "Axis Font", el("box plot viewer"), "normal normal 10px \"Roboto\""));
      await test.step("Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Date category mapping", async () => {
      await test.step("When user sets \"Category 1\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Category 1", el("box plot viewer"), "STARTED"));
      await test.step("Then \"Category 1\" property of box plot viewer should be \"STARTED\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "STARTED"));
      await test.step("When user sets \"Category 1 Map\" property of box plot viewer to \"month\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "month"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("And \"Category 1 Map\" property of box plot viewer should be \"month\"", () => propertyShouldBe(page, "Category 1 Map", el("box plot viewer"), "month"));
      await test.step("When user sets \"Category 1 Map\" property of box plot viewer to \"quarter\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "quarter"));
      await test.step("Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await test.step("When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await test.step("Then \"Category 1\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "RACE"));
      await test.step("And no errors should have been logged", () => noErrors(page));
      await test.step("When user sets \"Category 1\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Category 1", el("box plot viewer"), "SEX"));
    });
    await run.scenario("Custom tooltip", async () => {
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","RACE"],["Marker Size","10"],["Row Tooltip","AGE\nSEX\nWEIGHT"],["Show Tooltip","show custom tooltip"]]));
      await test.step("Then \"Row Tooltip\" property of box plot viewer should be \"AGE\\nSEX\\nWEIGHT\"", () => propertyShouldBe(page, "Row Tooltip", el("box plot viewer"), "AGE\\nSEX\\nWEIGHT"));
      await test.step("When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await test.step("Then the tooltip should show columns \"AGE, SEX, WEIGHT\"", () => tooltipColumns(page, "AGE, SEX, WEIGHT"));
      await test.step("When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await test.step("Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Category 1","SEX"]]));
      await test.step("And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await test.step("Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await test.step("And the tooltip should not show columns \"AGE, SEX, WEIGHT\"", () => tooltipNotColumns(page, "AGE, SEX, WEIGHT"));
      await test.step("When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
    });
    await run.scenario("Table switching resets Category 2", async () => {
      await test.step("Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await test.step("And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await test.step("When user sets \"Category 2\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 2", el("box plot viewer"), "RACE"));
      await test.step("Then \"Category 2\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Category 2", el("box plot viewer"), "RACE"));
      await test.step("When user sets \"Table\" property of box plot viewer to \"spgi-100\"", () => setProperty(page, "Table", el("box plot viewer"), "spgi-100"));
      await test.step("Then \"Category 2\" property of box plot viewer should not be \"RACE\"", () => propertyShouldNotBe(page, "Category 2", el("box plot viewer"), "RACE"));
      await test.step("And no errors should have been logged", () => noErrors(page));
      await test.step("When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","Average Mass"],["Category 1","Series"]]));
      await test.step("Then \"Table\" property of box plot viewer should be \"spgi-100\"", () => propertyShouldBe(page, "Table", el("box plot viewer"), "spgi-100"));
      await test.step("And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await test.step("When user sets \"Table\" property of box plot viewer to \"demog-1000\"", () => setProperty(page, "Table", el("box plot viewer"), "demog-1000"));
      await test.step("And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","AGE"],["Category 1","SEX"],["Category 2",""]]));
      await test.step("Then \"Table\" property of box plot viewer should be \"demog-1000\"", () => propertyShouldBe(page, "Table", el("box plot viewer"), "demog-1000"));
    });
    await run.scenario("Coloring keeps the render valid", async () => {
      await test.step("When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await test.step("Then \"Marker Color Column\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await test.step("And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await test.step("And no errors should have been logged", () => noErrors(page));
      await test.step("When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Double-click resets the view", async () => {
      await test.step("Given user listens for \"d4-boxplot-reset-view\" event on box plot viewer", () => listenFor(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await test.step("When user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
      await test.step("Then box plot viewer should show a narrowed value range", () => narrowedRange(page, el("box plot viewer")));
      await test.step("When user double-clicks on empty plot space of box plot viewer", () => doubleClickEmptySpace(page, el("box plot viewer")));
      await test.step("Then \"d4-boxplot-reset-view\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await test.step("And box plot viewer should show the full value range again", () => fullRange(page, el("box plot viewer")));
    });
    run.finish();
  });
});
