/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/ligand-rows.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.biostructure, biostructureviewer.viewer.ngl]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, isExpanded, pressKeyIn, selectIn, shouldBe, typeInto, uncheck, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {mouseOverRowIs, noneSelected, rowsRangeSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, clickAreaHolding, hoverArea, noBalloons, noErrors, propertyShouldBe, readingIs, readingReads, readingSame, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Which rows' ligands the Biostructure and NGL viewers show — selected, current and mouse-over rows", () => {
  const session = feature(test, "features/biostructure-viewer/ligand-rows.feature", import.meta.url);
  test("Which rows' ligands the Biostructure and NGL viewers show — selected, current and mouse-over rows", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.viewer.biostructure", "@realizes:biostructureviewer.viewer.ngl"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(24, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(25, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(26, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(27, "When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree")));
    await session.step(28, "Then the \"1bdq-obs-pred\" table view should open with 22 rows", () => tableViewOpened(page, "1bdq-obs-pred", 22));
    await session.step(29, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(30, "And user types \"Biostructure\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Biostructure", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(31, "And user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
    await session.step(32, "Then \"Data File\" input in Biostructure viewer should be visible", () => shouldBe(page, el("\"Data File\" input in Biostructure viewer"), "visible"));
    await session.step(33, "When user clicks on \"folder-tree\" icon in \"Data File\" input in Biostructure viewer", () => clickOn(page, el("\"folder-tree\" icon in \"Data File\" input in Biostructure viewer")));
    await session.step(34, "Then \"Select a file\" dialog should be visible", () => shouldBe(page, el("\"Select a file\" dialog"), "visible"));
    await session.step(35, "Given Files---App-Data tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("Files---App-Data tree node in \"Select a file\" dialog")));
    await session.step(36, "And Files---App-Data---BiostructureViewer tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node in \"Select a file\" dialog")));
    await session.step(37, "And Files---App-Data---BiostructureViewer---samples tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node in \"Select a file\" dialog")));
    await session.step(38, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node in \"Select a file\" dialog", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node in \"Select a file\" dialog")));
    await session.step(39, "And user clicks on OK button in \"Select a file\" dialog", () => clickOn(page, el("OK button in \"Select a file\" dialog")));
    await session.step(40, "Then \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
    await session.step(41, "And \"ligandColumnName\" property of Biostructure viewer should be \"molecule\"", () => propertyShouldBe(page, "ligandColumnName", el("Biostructure viewer"), "molecule"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Biostructure shows the ligands of the selected rows with Show Selected Rows Ligands", async () => {
      await session.step(46, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
      await session.step(47, "Given \"Behaviour\" category in context panel is expanded", () => isExpanded(page, el("\"Behaviour\" category in context panel")));
      await session.step(48, "When user unchecks \"Show Current Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Current Row Ligand\" property in context panel")));
      await session.step(49, "And user unchecks \"Show Mouse-Over Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Mouse-Over Row Ligand\" property in context panel")));
      await session.step(50, "Then the \"ligands shown\" reading of Biostructure viewer should be 0", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 0));
      await session.step(51, "When user takes a snapshot of Biostructure viewer", () => takeSnapshot(page, el("Biostructure viewer")));
      await session.step(52, "And user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")));
      await session.step(53, "And user clicks on the \"row header 4\" area of grid holding Shift", () => clickAreaHolding(page, "row header 4", el("grid"), "Shift"));
      await session.step(54, "Then rows 1 to 4 should be selected", () => rowsRangeSelected(page, 1, 4));
      await session.step(55, "And 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(56, "And the \"ligands shown\" reading of Biostructure viewer should be the same as before", () => readingSame(page, "ligands shown", el("Biostructure viewer")));
      await session.step(57, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
      await session.step(58, "And user checks \"Show Selected Rows Ligands\" property in context panel", () => check(page, el("\"Show Selected Rows Ligands\" property in context panel")));
      await session.step(59, "Then the \"ligand rows\" reading of Biostructure viewer should be \"1, 2, 3, 4\"", () => readingReads(page, "ligand rows", el("Biostructure viewer"), "1, 2, 3, 4"));
      await session.step(60, "And the \"ligands shown\" reading of Biostructure viewer should be 4", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 4));
      await session.step(61, "When user presses Escape in grid", () => pressKeyIn(page, "Escape", el("grid")));
      await session.step(62, "Then no rows should be selected", () => noneSelected(page));
      await session.step(63, "And the \"ligands shown\" reading of Biostructure viewer should be 0", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 0));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("NGL shows the same selected rows' ligands as Biostructure (GROK-17967)", async () => {
      await session.step(68, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(69, "And user types \"NGL\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "NGL", el("viewer gallery search in \"Add Viewer\" dialog")));
      await session.step(70, "And user clicks on first \"NGL\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"NGL\" button in \"Add Viewer\" dialog")));
      await session.step(71, "Then NGL viewer should be visible", () => shouldBe(page, el("NGL viewer"), "visible"));
      await session.step(72, "When user uploads \"../../BiostructureViewer/files/samples/1bdq.pdb\" through \"Open...\" link in NGL viewer", () => uploadThrough(page, "../../BiostructureViewer/files/samples/1bdq.pdb", el("\"Open...\" link in NGL viewer")));
      await session.step(73, "Then the \"structure loaded\" reading of NGL viewer should be \"true\"", () => readingReads(page, "structure loaded", el("NGL viewer"), "true"));
      await session.step(74, "When user clicks on settings icon of NGL viewer", () => clickOn(page, el("settings icon of NGL viewer")));
      await session.step(75, "Given \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(76, "When user selects \"molecule\" in \"Ligand\" property in context panel", () => selectIn(page, "molecule", el("\"Ligand\" property in context panel")));
      await session.step(77, "Then \"ligandColumnName\" property of NGL viewer should be \"molecule\"", () => propertyShouldBe(page, "ligandColumnName", el("NGL viewer"), "molecule"));
      await session.step(78, "Given \"Behaviour\" category in context panel is expanded", () => isExpanded(page, el("\"Behaviour\" category in context panel")));
      await session.step(79, "When user unchecks \"Show Current Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Current Row Ligand\" property in context panel")));
      await session.step(80, "And user unchecks \"Show Mouse-Over Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Mouse-Over Row Ligand\" property in context panel")));
      await session.step(81, "And user checks \"Show Selected Rows Ligands\" property in context panel", () => check(page, el("\"Show Selected Rows Ligands\" property in context panel")));
      await session.step(82, "Then the \"ligands shown\" reading of NGL viewer should be 0", () => readingIs(page, "ligands shown", el("NGL viewer"), 0));
      await session.step(83, "When user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")));
      await session.step(84, "And user clicks on the \"row header 4\" area of grid holding Shift", () => clickAreaHolding(page, "row header 4", el("grid"), "Shift"));
      await session.step(85, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(86, "And the \"ligand rows\" reading of NGL viewer should be \"1, 2, 3, 4\"", () => readingReads(page, "ligand rows", el("NGL viewer"), "1, 2, 3, 4"));
      await session.step(87, "And the \"ligands shown\" reading of NGL viewer should be 4", () => readingIs(page, "ligands shown", el("NGL viewer"), 4));
      await session.step(88, "And the \"ligand rows\" reading of Biostructure viewer should be \"1, 2, 3, 4\"", () => readingReads(page, "ligand rows", el("Biostructure viewer"), "1, 2, 3, 4"));
      await session.step(89, "When user presses Escape in grid", () => pressKeyIn(page, "Escape", el("grid")));
      await session.step(90, "Then no rows should be selected", () => noneSelected(page));
      await session.step(91, "And the \"ligands shown\" reading of NGL viewer should be 0", () => readingIs(page, "ligands shown", el("NGL viewer"), 0));
      await session.step(92, "And the \"ligands shown\" reading of Biostructure viewer should be 0", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 0));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
      await session.step(94, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("NGL's Show Current Row Ligand and Show Mouse-Over Row Ligand switch their ligand on and off", async () => {
      await session.step(97, "When user clicks on settings icon of NGL viewer", () => clickOn(page, el("settings icon of NGL viewer")));
      await session.step(98, "Given \"Behaviour\" category in context panel is expanded", () => isExpanded(page, el("\"Behaviour\" category in context panel")));
      await session.step(99, "When user checks \"Show Current Row Ligand\" property in context panel", () => check(page, el("\"Show Current Row Ligand\" property in context panel")));
      await session.step(100, "And user clicks on the \"cell 6 of molecule\" area of grid", () => clickArea(page, "cell 6 of molecule", el("grid")));
      await session.step(101, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(102, "And the \"ligand rows\" reading of NGL viewer should be \"6\"", () => readingReads(page, "ligand rows", el("NGL viewer"), "6"));
      await session.step(103, "And the \"ligands shown\" reading of NGL viewer should be 1", () => readingIs(page, "ligands shown", el("NGL viewer"), 1));
      await session.step(104, "When user clicks on settings icon of NGL viewer", () => clickOn(page, el("settings icon of NGL viewer")));
      await session.step(105, "And user unchecks \"Show Current Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Current Row Ligand\" property in context panel")));
      await session.step(106, "Then the \"ligands shown\" reading of NGL viewer should be 0", () => readingIs(page, "ligands shown", el("NGL viewer"), 0));
      await session.step(107, "When user takes a snapshot of NGL viewer", () => takeSnapshot(page, el("NGL viewer")));
      await session.step(108, "And user hovers over the \"cell 9 of molecule\" area of grid", () => hoverArea(page, "cell 9 of molecule", el("grid")));
      await session.step(109, "Then the mouse-over row of the table should be 9", () => mouseOverRowIs(page, 9));
      await session.step(110, "And the \"ligands shown\" reading of NGL viewer should be the same as before", () => readingSame(page, "ligands shown", el("NGL viewer")));
      await session.step(111, "When user clicks on settings icon of NGL viewer", () => clickOn(page, el("settings icon of NGL viewer")));
      await session.step(112, "And user checks \"Show Mouse-Over Row Ligand\" property in context panel", () => check(page, el("\"Show Mouse-Over Row Ligand\" property in context panel")));
      await session.step(113, "And user hovers over the \"cell 9 of molecule\" area of grid", () => hoverArea(page, "cell 9 of molecule", el("grid")));
      await session.step(114, "Then the mouse-over row of the table should be 9", () => mouseOverRowIs(page, 9));
      await session.step(115, "And the \"ligand rows\" reading of NGL viewer should be \"9\"", () => readingReads(page, "ligand rows", el("NGL viewer"), "9"));
      await session.step(116, "And the \"ligands shown\" reading of NGL viewer should be 1", () => readingIs(page, "ligands shown", el("NGL viewer"), 1));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
      await session.step(118, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
