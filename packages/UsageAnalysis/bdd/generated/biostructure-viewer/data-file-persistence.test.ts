/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/data-file-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.biostructure, views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, typeInto, uncheck, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, dialogCloses, noProjectOnServer, projectsOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, infoBalloonText, noBalloons, noErrors, pickFromContextMenu, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A structure loaded into an empty Biostructure viewer survives saving and reopening the project", () => {
  const session = feature(test, "features/biostructure-viewer/data-file-persistence.feature", import.meta.url);
  test("A structure loaded into an empty Biostructure viewer survives saving and reopening the project", {tag: ["@journey", "@serial", "@viewers", "@realizes:biostructureviewer.viewer.biostructure", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And simple mode is off", () => simpleModeOff(page));
    await session.step(23, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(24, "And no project named \"BsvDataFilePersistence\" is on the server", () => noProjectOnServer(page, "BsvDataFilePersistence"));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(27, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(28, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(29, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await run.scenario("The empty Biostructure viewer loads a structure chosen in its Data file input (GROK-16143)", async () => {
      await session.step(32, "When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree")));
      await session.step(33, "Then the \"1bdq-obs-pred\" table view should open with 22 rows", () => tableViewOpened(page, "1bdq-obs-pred", 22));
      await session.step(34, "And \"molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecule", "Molecule"));
      await session.step(35, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
      await session.step(36, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(37, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
      await session.step(38, "When user types \"Biostructure\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Biostructure", el("viewer gallery search in \"Add Viewer\" dialog")));
      await session.step(39, "And user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
      await session.step(40, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
      await session.step(41, "And \"Data File\" input in Biostructure viewer should be visible", () => shouldBe(page, el("\"Data File\" input in Biostructure viewer"), "visible"));
      await session.step(42, "And \"Reset Camera\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "absent"));
      await session.step(43, "When user clicks on \"folder-tree\" icon in \"Data File\" input in Biostructure viewer", () => clickOn(page, el("\"folder-tree\" icon in \"Data File\" input in Biostructure viewer")));
      await session.step(44, "Then \"Select a file\" dialog should be visible", () => shouldBe(page, el("\"Select a file\" dialog"), "visible"));
      await session.step(45, "Given Files---App-Data tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("Files---App-Data tree node in \"Select a file\" dialog")));
      await session.step(46, "And Files---App-Data---BiostructureViewer tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node in \"Select a file\" dialog")));
      await session.step(47, "And Files---App-Data---BiostructureViewer---samples tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node in \"Select a file\" dialog")));
      await session.step(48, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node in \"Select a file\" dialog", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node in \"Select a file\" dialog")));
      await session.step(49, "And user clicks on OK button in \"Select a file\" dialog", () => clickOn(page, el("OK button in \"Select a file\" dialog")));
      await session.step(50, "Then \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(51, "And \"Data File\" input in Biostructure viewer should be absent", () => shouldBe(page, el("\"Data File\" input in Biostructure viewer"), "absent"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Only the current row's ligand is shown, in Biostructure and in NGL (GROK-17967)", async () => {
      await session.step(57, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
      await session.step(58, "Given \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(59, "When user selects \"molecule\" in \"Ligand\" property in context panel", () => selectIn(page, "molecule", el("\"Ligand\" property in context panel")));
      await session.step(60, "Then \"ligandColumnName\" property of Biostructure viewer should be \"molecule\"", () => propertyShouldBe(page, "ligandColumnName", el("Biostructure viewer"), "molecule"));
      await session.step(61, "When user clicks on the \"cell 1 of molecule\" area of grid", () => clickArea(page, "cell 1 of molecule", el("grid")));
      await session.step(62, "Then the \"ligand rows\" reading of Biostructure viewer should be \"1\"", () => readingReads(page, "ligand rows", el("Biostructure viewer"), "1"));
      await session.step(63, "And the \"ligands shown\" reading of Biostructure viewer should be 1", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 1));
      await session.step(64, "When user clicks on the \"cell 5 of molecule\" area of grid", () => clickArea(page, "cell 5 of molecule", el("grid")));
      await session.step(65, "Then the \"ligand rows\" reading of Biostructure viewer should be \"5\"", () => readingReads(page, "ligand rows", el("Biostructure viewer"), "5"));
      await session.step(66, "And the \"ligands shown\" reading of Biostructure viewer should be 1", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 1));
      await session.step(67, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
      await session.step(68, "Given \"Behaviour\" category in context panel is expanded", () => isExpanded(page, el("\"Behaviour\" category in context panel")));
      await session.step(69, "When user unchecks \"Show Current Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Current Row Ligand\" property in context panel")));
      await session.step(70, "And user unchecks \"Show Mouse-Over Row Ligand\" property in context panel", () => uncheck(page, el("\"Show Mouse-Over Row Ligand\" property in context panel")));
      await session.step(71, "Then the \"ligands shown\" reading of Biostructure viewer should be 0", () => readingIs(page, "ligands shown", el("Biostructure viewer"), 0));
      await session.step(72, "When user checks \"Show Current Row Ligand\" property in context panel", () => check(page, el("\"Show Current Row Ligand\" property in context panel")));
      await session.step(73, "And user checks \"Show Mouse-Over Row Ligand\" property in context panel", () => check(page, el("\"Show Mouse-Over Row Ligand\" property in context panel")));
      await session.step(74, "Then the \"ligand rows\" reading of Biostructure viewer should be \"5\"", () => readingReads(page, "ligand rows", el("Biostructure viewer"), "5"));
      await session.step(75, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(76, "And user types \"NGL\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "NGL", el("viewer gallery search in \"Add Viewer\" dialog")));
      await session.step(77, "And user clicks on first \"NGL\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"NGL\" button in \"Add Viewer\" dialog")));
      await session.step(78, "Then NGL viewer should be visible", () => shouldBe(page, el("NGL viewer"), "visible"));
      await session.step(79, "When user uploads \"../../BiostructureViewer/files/samples/1bdq.pdb\" through \"Open...\" link in NGL viewer", () => uploadThrough(page, "../../BiostructureViewer/files/samples/1bdq.pdb", el("\"Open...\" link in NGL viewer")));
      await session.step(80, "Then the \"structure loaded\" reading of NGL viewer should be \"true\"", () => readingReads(page, "structure loaded", el("NGL viewer"), "true"));
      await session.step(81, "When user clicks on settings icon of NGL viewer", () => clickOn(page, el("settings icon of NGL viewer")));
      await session.step(82, "Given \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(83, "When user selects \"molecule\" in \"Ligand\" property in context panel", () => selectIn(page, "molecule", el("\"Ligand\" property in context panel")));
      await session.step(84, "Then \"ligandColumnName\" property of NGL viewer should be \"molecule\"", () => propertyShouldBe(page, "ligandColumnName", el("NGL viewer"), "molecule"));
      await session.step(85, "When user clicks on the \"cell 3 of molecule\" area of grid", () => clickArea(page, "cell 3 of molecule", el("grid")));
      await session.step(86, "Then the \"ligand rows\" reading of NGL viewer should be \"3\"", () => readingReads(page, "ligand rows", el("NGL viewer"), "3"));
      await session.step(87, "And the \"ligands shown\" reading of NGL viewer should be 1", () => readingIs(page, "ligands shown", el("NGL viewer"), 1));
      await session.step(88, "And the \"ligand rows\" reading of Biostructure viewer should be \"3\"", () => readingReads(page, "ligand rows", el("Biostructure viewer"), "3"));
      await session.step(89, "When user clicks on close icon of NGL viewer", () => clickOn(page, el("close icon of NGL viewer")));
      await session.step(90, "Then NGL viewer should be absent", () => shouldBe(page, el("NGL viewer"), "absent"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
      await session.step(92, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The structure comes back when the project is saved and reopened (GROK-17485)", async () => {
      await session.step(95, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(96, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(97, "When user enters \"BsvDataFilePersistence\" into Name text input in \"Save project\" dialog", () => enterInto(page, "BsvDataFilePersistence", el("Name text input in \"Save project\" dialog")));
      await session.step(98, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(99, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(100, "And an info balloon containing 'Project \"BsvDataFilePersistence\" uploaded' should have been shown", () => infoBalloonText(page, "Project \"BsvDataFilePersistence\" uploaded"));
      await session.step(101, "And \"Share BsvDataFilePersistence\" dialog should be visible", () => shouldBe(page, el("\"Share BsvDataFilePersistence\" dialog"), "visible"));
      await session.step(102, "When user clicks on CANCEL button in \"Share BsvDataFilePersistence\" dialog", () => clickOn(page, el("CANCEL button in \"Share BsvDataFilePersistence\" dialog")));
      await session.step(103, "Then the \"Share BsvDataFilePersistence\" dialog should close", () => dialogCloses(page, "Share BsvDataFilePersistence"));
      await session.step(104, "And 1 project named \"BsvDataFilePersistence\" should be on the server", () => projectsOnServer(page, 1, "BsvDataFilePersistence"));
      await session.step(105, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(106, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(107, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(108, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(109, "And user enters \"BsvDataFilePersistence\" into gallery search", () => enterInto(page, "BsvDataFilePersistence", el("gallery search")));
      await session.step(110, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(111, "And user double-clicks on BsvDataFilePersistence gallery card", () => doubleClickOn(page, el("BsvDataFilePersistence gallery card")));
      await session.step(112, "Then the \"1bdq-obs-pred\" table view should open with 22 rows", () => tableViewOpened(page, "1bdq-obs-pred", 22));
      await session.step(113, "And \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(114, "And \"Data File\" input in Biostructure viewer should be absent", () => shouldBe(page, el("\"Data File\" input in Biostructure viewer"), "absent"));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
      await session.step(116, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
