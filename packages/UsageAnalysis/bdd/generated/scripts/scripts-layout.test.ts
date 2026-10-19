/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cleanLayouts} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, currentViewType, entityHasLayout, saveScript, scriptOnServer, scriptsView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasNoArea, noBalloons, noErrors, painted, pickFromAreaContextMenu, pickFromContextMenu, readingIs, readingsDiffer, readingsEqual, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The layout of a script's result", () => {
  const session = feature(test, "features/scripts/scripts-layout.feature", import.meta.url);
  test("The layout of a script's result", {tag: ["@journey", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And a script \"BddScriptLayout{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptLayout{time}"), "//language: javascript\n//output: dataframe df\ndf = await grok.data.getDemoTable('cars.csv');"));
    await session.step(26, "And the layouts saved for the script are deleted at the end", () => cleanLayouts(page));
    await session.step(27, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("The Layout tab asks for a run first", async () => {
      await session.step(30, "When user types \"BddScriptLayout{time}\" into gallery search", () => typeInto(page, session.text("BddScriptLayout{time}"), el("gallery search")));
      await session.step(31, "And user double-clicks on \"BddScriptLayout{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BddScriptLayout{time}\" link in gallery"))));
      await session.step(32, "Then the \"BddScriptLayout{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptLayout{time}")));
      await session.step(33, "When user clicks on Layout tab", () => clickOn(page, el("Layout tab")));
      await session.step(34, "Then \"Run script to get data and edit layout.\" text should be visible", () => shouldBe(page, el("\"Run script to get data and edit layout.\" text"), "visible"));
      await session.step(35, "And grid viewer should be absent", () => shouldBe(page, el("grid viewer"), "absent"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running the script fills the layout with its result", async () => {
      await session.step(39, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(40, "Then grid viewer should be visible", () => shouldBe(page, el("grid viewer"), "visible"));
      await session.step(41, "And the \"rows\" reading of grid should be 30", () => readingIs(page, "rows", el("grid"), 30));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A viewer added from the toolbox joins the layout", async () => {
      await session.step(46, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")));
      await session.step(47, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(48, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A colouring and a hidden column join the layout", async () => {
      await session.step(53, "Then the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(54, "When user picks \"Color Coding > Linear\" from the context menu of the \"header length\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header length", el("grid")));
      await session.step(55, "Then the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(56, "When user picks \"Hide\" from the context menu of the \"header model\" area of grid", () => pickFromAreaContextMenu(page, "Hide", "header model", el("grid")));
      await session.step(57, "Then grid should not have a \"header model\" area", () => hasNoArea(page, el("grid"), "header model"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Save stores the layout with the script", async () => {
      await session.step(61, "When user saves the script", () => saveScript(page));
      await session.step(62, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(63, "And the script \"BddScriptLayout{time}\" on the server should have a layout", () => entityHasLayout(page, "script", session.text("BddScriptLayout{time}")));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Run... from the Scripts view opens the result with the saved layout", async () => {
      await session.step(67, "When user closes all views", () => closeAllViews(page));
      await session.step(68, "Given user opens the Scripts view", () => scriptsView(page));
      await session.step(69, "When user types \"BddScriptLayout{time}\" into gallery search", () => typeInto(page, session.text("BddScriptLayout{time}"), el("gallery search")));
      await session.step(70, "And user picks \"Run...\" from the context menu of \"BddScriptLayout{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptLayout{time}\" link in gallery"))));
      await session.step(71, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(72, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(73, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(74, "And grid should not have a \"header model\" area", () => hasNoArea(page, el("grid"), "header model"));
      await session.step(75, "And the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
