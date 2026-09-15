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
import {clickOn, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendLists, legendSlot, loadLayout, noBalloons, noErrors, propertiesShouldBe, readingIs, saveLayoutToServer, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot legend and lifecycle", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-legend-and-lifecycle.feature", import.meta.url);
  test("Trellis plot legend and lifecycle", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(22, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await run.scenario("A colour column gives the trellis a legend", async () => {
      await session.step(25, "Then legend of trellis plot viewer should be hidden", () => shouldBe(page, el("legend of trellis plot viewer"), "hidden"));
      await session.step(26, "When user sets \"colorColumnName\" inner property of trellis plot viewer to \"SEX\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "SEX"));
      await session.step(27, "And user sets \"Legend Visibility\" property of trellis plot viewer to \"Always\"", () => setProperty(page, "Legend Visibility", el("trellis plot viewer"), "Always"));
      await session.step(28, "Then legend of trellis plot viewer should be visible", () => shouldBe(page, el("legend of trellis plot viewer"), "visible"));
      await session.step(29, "And the legend of trellis plot viewer should list 2 items", () => legendLists(page, el("trellis plot viewer"), 2));
      await session.step(30, "When user sets \"colorColumnName\" inner property of trellis plot viewer to \"RACE\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "RACE"));
      await session.step(31, "Then the legend of trellis plot viewer should list 4 items", () => legendLists(page, el("trellis plot viewer"), 4));
      await session.step(32, "And legend of trellis plot viewer should contain text \"Caucasian\"", () => shouldContainText(page, el("legend of trellis plot viewer"), "Caucasian"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The legend takes each of the four slots", async () => {
      await session.step(36, "When user sets \"Legend Position\" property of trellis plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Left"));
      await session.step(37, "Then the legend of trellis plot viewer should be in the \"left\" slot", () => legendSlot(page, el("trellis plot viewer"), "left"));
      await session.step(38, "When user sets \"Legend Position\" property of trellis plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Right"));
      await session.step(39, "Then the legend of trellis plot viewer should be in the \"right\" slot", () => legendSlot(page, el("trellis plot viewer"), "right"));
      await session.step(40, "When user sets \"Legend Position\" property of trellis plot viewer to \"Top\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Top"));
      await session.step(41, "Then the legend of trellis plot viewer should be in the \"top\" slot", () => legendSlot(page, el("trellis plot viewer"), "top"));
      await session.step(42, "When user sets \"Legend Position\" property of trellis plot viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("trellis plot viewer"), "Bottom"));
      await session.step(43, "Then the legend of trellis plot viewer should be in the \"bottom\" slot", () => legendSlot(page, el("trellis plot viewer"), "bottom"));
      await session.step(44, "And the legend of trellis plot viewer should list 4 items", () => legendLists(page, el("trellis plot viewer"), 4));
      await session.step(45, "When user sets \"Legend Visibility\" property of trellis plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("trellis plot viewer"), "Never"));
      await session.step(46, "Then legend of trellis plot viewer should be hidden", () => shouldBe(page, el("legend of trellis plot viewer"), "hidden"));
      await session.step(47, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Legend Visibility","Always"],["Legend Position","Auto"]]));
      await session.step(50, "Then legend of trellis plot viewer should be visible", () => shouldBe(page, el("legend of trellis plot viewer"), "visible"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Undo brings the closed viewer back and redo closes it again", async () => {
      await session.step(54, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(55, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(56, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(57, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(58, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(59, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(60, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(61, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(62, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(63, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(64, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout round-trip through the server restores the split and the type", async () => {
      await session.step(69, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]));
      await session.step(73, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(74, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(75, "And user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(76, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(77, "When user loads the saved layout", () => loadLayout(page));
      await session.step(78, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(79, "And properties of trellis plot viewer should be:", () => propertiesShouldBe(page, el("trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEX"],["Viewer Type","Bar chart"]]));
      await session.step(83, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(84, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip restores the split, the type and the row source", async () => {
      await session.step(88, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"],["Row Source","All"]]));
      await session.step(93, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(94, "When user saves the current view as project \"bdd-trellis-plot\"", () => saveAsProject(page, "bdd-trellis-plot"));
      await session.step(95, "And user closes all views", () => closeAllViews(page));
      await session.step(96, "And user opens the \"bdd-trellis-plot\" project", () => openProject(page, "bdd-trellis-plot"));
      await session.step(97, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(98, "And properties of trellis plot viewer should be:", () => propertiesShouldBe(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"],["Row Source","All"]]));
      await session.step(103, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(104, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(105, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
