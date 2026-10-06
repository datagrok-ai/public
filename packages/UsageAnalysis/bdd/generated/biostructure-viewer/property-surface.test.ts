/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/property-surface.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.biostructure]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, isExpanded, pressKey, selectIn, shouldBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, contextPanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noBalloons, noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Biostructure viewer settings — Binding Site Whole Residues and the Controls category", () => {
  const session = feature(test, "features/biostructure-viewer/property-surface.feature", import.meta.url);
  test("Biostructure viewer settings — Binding Site Whole Residues and the Controls category", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.viewer.biostructure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(28, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(29, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(30, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(31, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(32, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(33, "When user clicks on the \"cell 1 of pdb\" area of grid", () => clickArea(page, "cell 1 of pdb", el("grid")));
    await session.step(34, "Then \"PDB Information\" section in context panel should be visible", () => shouldBe(page, el("\"PDB Information\" section in context panel"), "visible"));
    await session.step(35, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
    await session.step(36, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(37, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(38, "When user types \"Biostructure\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Biostructure", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(39, "And user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
    await session.step(40, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(41, "And Biostructure viewer should be visible", () => shouldBe(page, el("Biostructure viewer"), "visible"));
    await session.step(42, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
    await session.step(43, "And user selects \"pdb\" in \"Biostructure Id\" property in context panel", () => selectIn(page, "pdb", el("\"Biostructure Id\" property in context panel")));
    await session.step(44, "Then \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Binding Site Whole Residues switches off and on under Show Binding Site", async () => {
      await session.step(49, "Given \"Binding Site\" category in context panel is expanded", () => isExpanded(page, el("\"Binding Site\" category in context panel")));
      await session.step(50, "Then \"Show Binding Site\" property in context panel should be unchecked", () => shouldBe(page, el("\"Show Binding Site\" property in context panel"), "unchecked"));
      await session.step(51, "And \"Binding Site Whole Residues\" property in context panel should be checked", () => shouldBe(page, el("\"Binding Site Whole Residues\" property in context panel"), "checked"));
      await session.step(52, "When user checks \"Show Binding Site\" property in context panel", () => check(page, el("\"Show Binding Site\" property in context panel")));
      await session.step(53, "Then \"showBindingSite\" property of Biostructure viewer should be \"true\"", () => propertyShouldBe(page, "showBindingSite", el("Biostructure viewer"), "true"));
      await session.step(54, "And \"Binding site\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Binding site\" overlay button in Biostructure viewer"), "selected"));
      await session.step(55, "When user clicks on \"Binding site\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Binding site\" overlay button in Biostructure viewer")));
      await session.step(56, "Then \"Show side chains\" checkbox should be visible", () => shouldBe(page, el("\"Show side chains\" checkbox"), "visible"));
      await session.step(57, "And \"Show side chains\" checkbox should be checked", () => shouldBe(page, el("\"Show side chains\" checkbox"), "checked"));
      await session.step(58, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(59, "Then \"Show side chains\" checkbox should be hidden", () => shouldBe(page, el("\"Show side chains\" checkbox"), "hidden"));
      await session.step(60, "When user unchecks \"Binding Site Whole Residues\" property in context panel", () => uncheck(page, el("\"Binding Site Whole Residues\" property in context panel")));
      await session.step(61, "Then \"bindingSiteWholeResidues\" property of Biostructure viewer should be \"false\"", () => propertyShouldBe(page, "bindingSiteWholeResidues", el("Biostructure viewer"), "false"));
      await session.step(62, "And \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(63, "When user checks \"Binding Site Whole Residues\" property in context panel", () => check(page, el("\"Binding Site Whole Residues\" property in context panel")));
      await session.step(64, "Then \"bindingSiteWholeResidues\" property of Biostructure viewer should be \"true\"", () => propertyShouldBe(page, "bindingSiteWholeResidues", el("Biostructure viewer"), "true"));
      await session.step(65, "When user unchecks \"Show Binding Site\" property in context panel", () => uncheck(page, el("\"Show Binding Site\" property in context panel")));
      await session.step(66, "Then \"showBindingSite\" property of Biostructure viewer should be \"false\"", () => propertyShouldBe(page, "showBindingSite", el("Biostructure viewer"), "false"));
      await session.step(67, "When user clicks on \"Binding site\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Binding site\" overlay button in Biostructure viewer")));
      await session.step(68, "Then \"Show side chains\" checkbox should be visible", () => shouldBe(page, el("\"Show side chains\" checkbox"), "visible"));
      await session.step(69, "And \"Show side chains\" checkbox should be unchecked", () => shouldBe(page, el("\"Show side chains\" checkbox"), "unchecked"));
      await session.step(70, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(71, "Then \"Show side chains\" checkbox should be hidden", () => shouldBe(page, el("\"Show side chains\" checkbox"), "hidden"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Controls category holds Show Welcome Toast and Show Import Controls, both off", async () => {
      await session.step(76, "Given \"Controls\" category in context panel is expanded", () => isExpanded(page, el("\"Controls\" category in context panel")));
      await session.step(77, "Then \"Show Welcome Toast\" property in context panel should be unchecked", () => shouldBe(page, el("\"Show Welcome Toast\" property in context panel"), "unchecked"));
      await session.step(78, "And \"Show Import Controls\" property in context panel should be unchecked", () => shouldBe(page, el("\"Show Import Controls\" property in context panel"), "unchecked"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
