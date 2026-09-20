/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-legend-and-lifecycle.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
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
import {clickOn, pressKey, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendCross, clickLegendItem, clickLegendItemHolding, legendLists, legendSlot, loadLayout, noBalloons, noErrors, propertiesShouldBe, readingIs, readingReads, saveLayoutToServer, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, innerPropertyShouldBe, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot legend and lifecycle", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-legend-and-lifecycle.feature", import.meta.url);
  test("Trellis plot legend and lifecycle", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(28, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await run.scenario("A colour column gives the trellis a legend", async () => {
      await session.step(31, "Then legend of trellis plot viewer should be hidden", () => shouldBe(page, el("legend of trellis plot viewer"), "hidden"));
      await session.step(32, "When user sets \"colorColumnName\" inner property of trellis plot viewer to \"SEX\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "SEX"));
      await session.step(33, "And user sets \"Legend Visibility\" property of trellis plot viewer to \"Always\"", () => setProperty(page, "Legend Visibility", el("trellis plot viewer"), "Always"));
      await session.step(34, "Then legend of trellis plot viewer should be visible", () => shouldBe(page, el("legend of trellis plot viewer"), "visible"));
      await session.step(35, "And the legend of trellis plot viewer should list 2 items", () => legendLists(page, el("trellis plot viewer"), 2));
      await session.step(36, "When user sets \"colorColumnName\" inner property of trellis plot viewer to \"RACE\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "RACE"));
      await session.step(37, "Then the legend of trellis plot viewer should list 4 items", () => legendLists(page, el("trellis plot viewer"), 4));
      await session.step(38, "And legend of trellis plot viewer should contain text \"Caucasian\"", () => shouldContainText(page, el("legend of trellis plot viewer"), "Caucasian"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The legend takes each of the four slots", async () => {
      await session.step(42, "When user sets \"Legend Position\" property of trellis plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Left"));
      await session.step(43, "Then the legend of trellis plot viewer should be in the \"left\" slot", () => legendSlot(page, el("trellis plot viewer"), "left"));
      await session.step(44, "When user sets \"Legend Position\" property of trellis plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Right"));
      await session.step(45, "Then the legend of trellis plot viewer should be in the \"right\" slot", () => legendSlot(page, el("trellis plot viewer"), "right"));
      await session.step(46, "When user sets \"Legend Position\" property of trellis plot viewer to \"Top\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Top"));
      await session.step(47, "Then the legend of trellis plot viewer should be in the \"top\" slot", () => legendSlot(page, el("trellis plot viewer"), "top"));
      await session.step(48, "When user sets \"Legend Position\" property of trellis plot viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Bottom"));
      await session.step(49, "Then the legend of trellis plot viewer should be in the \"bottom\" slot", () => legendSlot(page, el("trellis plot viewer"), "bottom"));
      await session.step(50, "And the legend of trellis plot viewer should list 4 items", () => legendLists(page, el("trellis plot viewer"), 4));
      await session.step(51, "When user sets \"Legend Visibility\" property of trellis plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("trellis plot viewer"), "Never"));
      await session.step(52, "Then legend of trellis plot viewer should be hidden", () => shouldBe(page, el("legend of trellis plot viewer"), "hidden"));
      await session.step(53, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Legend Visibility","Always"],["Legend Position","Auto"]]));
      await session.step(56, "Then legend of trellis plot viewer should be visible", () => shouldBe(page, el("legend of trellis plot viewer"), "visible"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A box plot inside keeps Show All Categories through the legend's clicks", async () => {
      await session.step(60, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Box plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Box plot"));
      await session.step(61, "And user sets \"showAllCategories\" inner property of trellis plot viewer to \"true\"", () => setInnerProperty(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(62, "Then \"showAllCategories\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(63, "And legend of trellis plot viewer should be visible", () => shouldBe(page, el("legend of trellis plot viewer"), "visible"));
      await session.step(64, "And the legend of trellis plot viewer should list 6 items", () => legendLists(page, el("trellis plot viewer"), 6));
      await session.step(65, "And \"AS\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"AS\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(66, "And \"Indigestion\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"Indigestion\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(67, "And \"PsA\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(68, "When user clicks on \"AS\" item in the legend of trellis plot viewer", () => clickLegendItem(page, "AS", el("trellis plot viewer")));
      await session.step(69, "Then \"AS\" legend item in legend of trellis plot viewer should be selected", () => shouldBe(page, el("\"AS\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(70, "And \"Indigestion\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"Indigestion\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(71, "And \"PsA\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(72, "And \"showAllCategories\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(73, "When user clicks on \"Indigestion\" item in the legend of trellis plot viewer holding Control", () => clickLegendItemHolding(page, "Indigestion", el("trellis plot viewer"), "Control"));
      await session.step(74, "Then \"AS\" legend item in legend of trellis plot viewer should be selected", () => shouldBe(page, el("\"AS\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(75, "And \"Indigestion\" legend item in legend of trellis plot viewer should be selected", () => shouldBe(page, el("\"Indigestion\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(76, "And \"PsA\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(77, "And \"showAllCategories\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(78, "When user clicks on the cross of \"AS\" item in the legend of trellis plot viewer", () => clickLegendCross(page, "AS", el("trellis plot viewer")));
      await session.step(79, "Then \"AS\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"AS\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(80, "And \"Indigestion\" legend item in legend of trellis plot viewer should be selected", () => shouldBe(page, el("\"Indigestion\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(81, "And \"PsA\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(82, "And \"showAllCategories\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(83, "When user clicks on \"PsA\" item in the legend of trellis plot viewer", () => clickLegendItem(page, "PsA", el("trellis plot viewer")));
      await session.step(84, "Then \"PsA\" legend item in legend of trellis plot viewer should be selected", () => shouldBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(85, "And \"Indigestion\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"Indigestion\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(86, "And \"AS\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"AS\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(87, "And \"showAllCategories\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(88, "When user clicks on the cross of \"PsA\" item in the legend of trellis plot viewer", () => clickLegendCross(page, "PsA", el("trellis plot viewer")));
      await session.step(89, "Then \"PsA\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(90, "And \"Indigestion\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"Indigestion\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(91, "And \"AS\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"AS\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(92, "And \"PsA\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"PsA\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(93, "And \"UC\" legend item in legend of trellis plot viewer should not be selected", () => shouldNotBe(page, el("\"UC\" legend item in legend of trellis plot viewer"), "selected"));
      await session.step(94, "And the legend of trellis plot viewer should list 6 items", () => legendLists(page, el("trellis plot viewer"), 6));
      await session.step(95, "And \"showAllCategories\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showAllCategories", el("trellis plot viewer"), "true"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
      await session.step(97, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(98, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("Undo brings the closed viewer back and redo closes it again", async () => {
      await session.step(101, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(102, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(103, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(104, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(105, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(106, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(107, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(108, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(109, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(110, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(111, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(112, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout round-trip through the server restores the split and the type", async () => {
      await session.step(116, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]));
      await session.step(120, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(121, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(122, "And user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(123, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(124, "When user loads the saved layout", () => loadLayout(page));
      await session.step(125, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(126, "And properties of trellis plot viewer should be:", () => propertiesShouldBe(page, el("trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]));
      await session.step(130, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(131, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip restores the split, the type and the row source", async () => {
      await session.step(135, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"],["Row Source","All"]]));
      await session.step(140, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(141, "When user saves the current view as project \"bdd-trellis-plot\"", () => saveAsProject(page, "bdd-trellis-plot"));
      await session.step(142, "And user closes all views", () => closeAllViews(page));
      await session.step(143, "And user opens the \"bdd-trellis-plot\" project", () => openProject(page, "bdd-trellis-plot"));
      await session.step(144, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(145, "And properties of trellis plot viewer should be:", () => propertiesShouldBe(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"],["Row Source","All"]]));
      await session.step(150, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(151, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(152, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(153, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
