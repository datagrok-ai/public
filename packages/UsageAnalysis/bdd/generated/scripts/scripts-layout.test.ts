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
import {clickOn, doubleClickOn, enterInto, pressKey, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, currentViewType, dialogCloses, entityHasLayout, noProjectOnServer, openProjectWithTable, projectsOnServer, saveScript, scriptOnServer, scriptsView, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickAreaHolding, hasNoArea, noBalloons, noErrors, painted, pickFromAreaContextMenu, pickFromContextMenu, readingIs, readingsDiffer, readingsEqual, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The layout of a script's result", () => {
  const session = feature(test, "features/scripts/scripts-layout.feature", import.meta.url);
  test("The layout of a script's result", {tag: ["@journey", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And a script \"BddScriptLayout{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptLayout{time}"), "//language: javascript\n//output: dataframe df\ndf = await grok.data.getDemoTable('cars.csv');"));
    await session.step(31, "And the layouts saved for the script are deleted at the end", () => cleanLayouts(page));
    await session.step(32, "And no project named \"BDD-Script-layout-proj-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Script-layout-proj-{time}")));
    await session.step(33, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("The Layout tab asks for a run first", async () => {
      await session.step(36, "When user types \"BddScriptLayout{time}\" into gallery search", () => typeInto(page, session.text("BddScriptLayout{time}"), el("gallery search")));
      await session.step(37, "And user double-clicks on \"BddScriptLayout{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BddScriptLayout{time}\" link in gallery"))));
      await session.step(38, "Then the \"BddScriptLayout{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptLayout{time}")));
      await session.step(39, "When user clicks on Layout tab", () => clickOn(page, el("Layout tab")));
      await session.step(40, "Then \"Run script to get data and edit layout.\" text should be visible", () => shouldBe(page, el("\"Run script to get data and edit layout.\" text"), "visible"));
      await session.step(41, "And grid viewer should be absent", () => shouldBe(page, el("grid viewer"), "absent"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running the script fills the layout with its result", async () => {
      await session.step(45, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(46, "Then grid viewer should be visible", () => shouldBe(page, el("grid viewer"), "visible"));
      await session.step(47, "And the \"rows\" reading of grid should be 30", () => readingIs(page, "rows", el("grid"), 30));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A viewer added from the toolbox joins the layout", async () => {
      await session.step(52, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")));
      await session.step(53, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(54, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A colouring and a hidden column join the layout", async () => {
      await session.step(59, "Then the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(60, "When user picks \"Color Coding > Linear\" from the context menu of the \"header length\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header length", el("grid")));
      await session.step(61, "Then the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(62, "When user picks \"Hide\" from the context menu of the \"header model\" area of grid", () => pickFromAreaContextMenu(page, "Hide", "header model", el("grid")));
      await session.step(63, "Then grid should not have a \"header model\" area", () => hasNoArea(page, el("grid"), "header model"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Save stores the layout with the script", async () => {
      await session.step(67, "When user saves the script", () => saveScript(page));
      await session.step(68, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(69, "And the script \"BddScriptLayout{time}\" on the server should have a layout", () => entityHasLayout(page, "script", session.text("BddScriptLayout{time}")));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Run... from the Scripts view opens the result with the saved layout", async () => {
      await session.step(73, "When user closes all views", () => closeAllViews(page));
      await session.step(74, "Given user opens the Scripts view", () => scriptsView(page));
      await session.step(75, "When user types \"BddScriptLayout{time}\" into gallery search", () => typeInto(page, session.text("BddScriptLayout{time}"), el("gallery search")));
      await session.step(76, "And user picks \"Run...\" from the context menu of \"BddScriptLayout{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptLayout{time}\" link in gallery"))));
      await session.step(77, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(78, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(79, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(80, "And grid should not have a \"header model\" area", () => hasNoArea(page, el("grid"), "header model"));
      await session.step(81, "And the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The result with one more viewer is saved as a project and opens with its layout", async () => {
      await session.step(86, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(87, "When user clicks on histogram icon in toolbox", () => clickOn(page, el("histogram icon in toolbox")));
      await session.step(88, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(89, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(90, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(91, "When user enters \"BDD-Script-layout-proj-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Script-layout-proj-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(92, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(93, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(94, "And 1 project named \"BDD-Script-layout-proj-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Script-layout-proj-{time}")));
      await session.step(95, "And \"Share BDD-Script-layout-proj-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Script-layout-proj-{time}\" dialog")), "visible"));
      await session.step(96, "When user clicks on CANCEL button in \"Share BDD-Script-layout-proj-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-Script-layout-proj-{time}\" dialog"))));
      await session.step(97, "Then the \"Share BDD-Script-layout-proj-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Script-layout-proj-{time}")));
      await session.step(98, "When user closes all views", () => closeAllViews(page));
      await session.step(99, "And user opens the \"BDD-Script-layout-proj-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Script-layout-proj-{time}")));
      await session.step(100, "Then the table should have 30 rows", () => rowCount(page, 30));
      await session.step(101, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(102, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(103, "And grid should not have a \"header model\" area", () => hasNoArea(page, el("grid"), "header model"));
      await session.step(104, "And the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
      await session.step(106, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("REFRESH of the reopened result runs the script again and keeps its layout", async () => {
      await session.step(110, "When user clicks on the \"row header 1\" area of grid holding Shift", () => clickAreaHolding(page, "row header 1", el("grid"), "Shift"));
      await session.step(111, "And user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(112, "Then the table should have 29 rows", () => rowCount(page, 29));
      await session.step(113, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(114, "When user clicks on REFRESH button in toolbox", () => clickOn(page, el("REFRESH button in toolbox")));
      await session.step(115, "Then the table should have 30 rows", () => rowCount(page, 30));
      await session.step(116, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(117, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(118, "And grid should not have a \"header model\" area", () => hasNoArea(page, el("grid"), "header model"));
      await session.step(119, "And the \"color of cell 1 of length\" and \"color of cell 2 of length\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of length", "color of cell 2 of length", el("grid")));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
      await session.step(121, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
