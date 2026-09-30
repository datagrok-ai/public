/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/data-file-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.biostructure, views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, dialogCloses, noProjectOnServer, projectsOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A structure loaded into an empty Biostructure viewer survives saving and reopening the project", () => {
  const session = feature(test, "features/biostructure-viewer/data-file-persistence.feature", import.meta.url);
  test("A structure loaded into an empty Biostructure viewer survives saving and reopening the project", {tag: ["@journey", "@serial", "@viewers", "@realizes:biostructureviewer.viewer.biostructure", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(25, "And no project named \"BsvDataFilePersistence\" is on the server", () => noProjectOnServer(page, "BsvDataFilePersistence"));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(28, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(29, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(30, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await run.scenario("The empty Biostructure viewer loads a structure chosen in its Data file input (GROK-16143)", async () => {
      await session.step(33, "When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq-obs-pred.sdf tree node inside browse tree")));
      await session.step(34, "Then the \"1bdq-obs-pred\" table view should open with 22 rows", () => tableViewOpened(page, "1bdq-obs-pred", 22));
      await session.step(35, "And \"molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecule", "Molecule"));
      await session.step(36, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
      await session.step(37, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
      await session.step(38, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
      await session.step(39, "When user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
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
    await run.scenario("The structure comes back when the project is saved and reopened (GROK-17485)", async () => {
      await session.step(56, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(57, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(58, "When user enters \"BsvDataFilePersistence\" into Name text input in \"Save project\" dialog", () => enterInto(page, "BsvDataFilePersistence", el("Name text input in \"Save project\" dialog")));
      await session.step(59, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(60, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(61, "And an info balloon containing 'Project \"BsvDataFilePersistence\" uploaded' should have been shown", () => infoBalloonText(page, "Project \"BsvDataFilePersistence\" uploaded"));
      await session.step(62, "And \"Share BsvDataFilePersistence\" dialog should be visible", () => shouldBe(page, el("\"Share BsvDataFilePersistence\" dialog"), "visible"));
      await session.step(63, "When user clicks on CANCEL button in \"Share BsvDataFilePersistence\" dialog", () => clickOn(page, el("CANCEL button in \"Share BsvDataFilePersistence\" dialog")));
      await session.step(64, "Then the \"Share BsvDataFilePersistence\" dialog should close", () => dialogCloses(page, "Share BsvDataFilePersistence"));
      await session.step(65, "And 1 project named \"BsvDataFilePersistence\" should be on the server", () => projectsOnServer(page, 1, "BsvDataFilePersistence"));
      await session.step(66, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(67, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(68, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(69, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(70, "And user enters \"BsvDataFilePersistence\" into gallery search", () => enterInto(page, "BsvDataFilePersistence", el("gallery search")));
      await session.step(71, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(72, "And user double-clicks on BsvDataFilePersistence gallery card", () => doubleClickOn(page, el("BsvDataFilePersistence gallery card")));
      await session.step(73, "Then the \"1bdq-obs-pred\" table view should open with 22 rows", () => tableViewOpened(page, "1bdq-obs-pred", 22));
      await session.step(74, "And \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(75, "And \"Data File\" input in Biostructure viewer should be absent", () => shouldBe(page, el("\"Data File\" input in Biostructure viewer"), "absent"));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
