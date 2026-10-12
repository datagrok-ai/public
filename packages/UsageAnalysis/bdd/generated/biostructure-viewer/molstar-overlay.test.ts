/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/molstar-overlay.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.biostructure]
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
import {check, clickOn, doubleClickOn, enterInto, isExpanded, pressKey, shouldBe, shouldNotBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Mol* viewport overlay buttons of the Biostructure viewer", () => {
  const session = feature(test, "features/biostructure-viewer/molstar-overlay.feature", import.meta.url);
  test("Mol* viewport overlay buttons of the Biostructure viewer", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.viewer.biostructure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(32, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(33, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(34, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(35, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(36, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(37, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
    await session.step(38, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(39, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(40, "When user types \"Biostructure\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Biostructure", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(41, "And user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
    await session.step(42, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(43, "And Biostructure viewer should be visible", () => shouldBe(page, el("Biostructure viewer"), "visible"));
    await session.step(44, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
    await session.step(45, "Then \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Screenshot / State Snapshot opens its panel and closes it again", async () => {
      await session.step(50, "Then \"Auto-crop\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Auto-crop\" button in Biostructure viewer"), "absent"));
      await session.step(51, "And \"Screenshot / State Snapshot\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Screenshot / State Snapshot\" overlay button in Biostructure viewer"), "selected"));
      await session.step(52, "When user clicks on \"Screenshot / State Snapshot\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Screenshot / State Snapshot\" overlay button in Biostructure viewer")));
      await session.step(53, "Then \"Auto-crop\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Auto-crop\" button in Biostructure viewer"), "visible"));
      await session.step(54, "And \"Screenshot / State Snapshot\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Screenshot / State Snapshot\" overlay button in Biostructure viewer"), "selected"));
      await session.step(55, "When user clicks on \"Screenshot / State Snapshot\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Screenshot / State Snapshot\" overlay button in Biostructure viewer")));
      await session.step(56, "Then \"Auto-crop\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Auto-crop\" button in Biostructure viewer"), "absent"));
      await session.step(57, "And \"Screenshot / State Snapshot\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Screenshot / State Snapshot\" overlay button in Biostructure viewer"), "selected"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toggle Controls Panel shows and hides the structure controls", async () => {
      await session.step(61, "Then \"Assembly\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "absent"));
      await session.step(62, "And \"Toggle Controls Panel\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer"), "selected"));
      await session.step(63, "When user clicks on \"Toggle Controls Panel\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer")));
      await session.step(64, "Then \"Assembly\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "visible"));
      await session.step(65, "And \"Toggle Controls Panel\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer"), "selected"));
      await session.step(66, "And the \"controls shown\" reading of Biostructure viewer should be \"true\"", () => readingReads(page, "controls shown", el("Biostructure viewer"), "true"));
      await session.step(67, "When user clicks on \"Toggle Controls Panel\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer")));
      await session.step(68, "Then \"Assembly\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "absent"));
      await session.step(69, "And \"Toggle Controls Panel\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer"), "selected"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Layout Show Controls in the settings shows and hides the same controls", async () => {
      await session.step(73, "Given \"Layout\" category in context panel is expanded", () => isExpanded(page, el("\"Layout\" category in context panel")));
      await session.step(74, "When user checks \"Layout Show Controls\" property in context panel", () => check(page, el("\"Layout Show Controls\" property in context panel")));
      await session.step(75, "Then \"Assembly\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "visible"));
      await session.step(76, "And \"Toggle Controls Panel\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer"), "selected"));
      await session.step(77, "When user unchecks \"Layout Show Controls\" property in context panel", () => uncheck(page, el("\"Layout Show Controls\" property in context panel")));
      await session.step(78, "Then \"Assembly\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "absent"));
      await session.step(79, "And \"Toggle Controls Panel\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Toggle Controls Panel\" overlay button in Biostructure viewer"), "selected"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toggle Selection Mode shows and hides the selection toolbar", async () => {
      await session.step(83, "Then \"Turn selection mode off\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Turn selection mode off\" button in Biostructure viewer"), "absent"));
      await session.step(84, "When user clicks on \"Toggle Selection Mode\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Toggle Selection Mode\" overlay button in Biostructure viewer")));
      await session.step(85, "Then \"Turn selection mode off\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Turn selection mode off\" button in Biostructure viewer"), "visible"));
      await session.step(86, "And \"Toggle Selection Mode\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Toggle Selection Mode\" overlay button in Biostructure viewer"), "selected"));
      await session.step(87, "When user clicks on \"Toggle Selection Mode\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Toggle Selection Mode\" overlay button in Biostructure viewer")));
      await session.step(88, "Then \"Turn selection mode off\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Turn selection mode off\" button in Biostructure viewer"), "absent"));
      await session.step(89, "And \"Toggle Selection Mode\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Toggle Selection Mode\" overlay button in Biostructure viewer"), "selected"));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Settings / Controls Info opens the Mol* settings panel and closes it again", async () => {
      await session.step(93, "Then \"Mouse Controls\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Mouse Controls\" button in Biostructure viewer"), "absent"));
      await session.step(94, "When user clicks on \"Settings / Controls Info\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Settings / Controls Info\" overlay button in Biostructure viewer")));
      await session.step(95, "Then \"Mouse Controls\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Mouse Controls\" button in Biostructure viewer"), "visible"));
      await session.step(96, "And \"Settings / Controls Info\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Settings / Controls Info\" overlay button in Biostructure viewer"), "selected"));
      await session.step(97, "When user clicks on \"Settings / Controls Info\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Settings / Controls Info\" overlay button in Biostructure viewer")));
      await session.step(98, "Then \"Mouse Controls\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Mouse Controls\" button in Biostructure viewer"), "absent"));
      await session.step(99, "And \"Settings / Controls Info\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Settings / Controls Info\" overlay button in Biostructure viewer"), "selected"));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Binding site button and the Binding Site settings follow each other", async () => {
      await session.step(103, "Then \"Binding site\" overlay button in Biostructure viewer should be enabled", () => shouldBe(page, el("\"Binding site\" overlay button in Biostructure viewer"), "enabled"));
      await session.step(104, "And \"Binding site\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Binding site\" overlay button in Biostructure viewer"), "selected"));
      await session.step(105, "And \"Binding Site\" text in dialog should be absent", () => shouldBe(page, el("\"Binding Site\" text in dialog"), "absent"));
      await session.step(106, "When user clicks on \"Binding site\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Binding site\" overlay button in Biostructure viewer")));
      await session.step(107, "Then \"Binding Site\" text in dialog should be visible", () => shouldBe(page, el("\"Binding Site\" text in dialog"), "visible"));
      await session.step(108, "And \"Show side chains\" checkbox should be visible", () => shouldBe(page, el("\"Show side chains\" checkbox"), "visible"));
      await session.step(109, "And \"Show side chains\" checkbox should be unchecked", () => shouldBe(page, el("\"Show side chains\" checkbox"), "unchecked"));
      await session.step(110, "And \"5.0 Å\" text should be visible", () => shouldBe(page, el("\"5.0 Å\" text"), "visible"));
      await session.step(111, "When user checks \"Show side chains\" checkbox", () => check(page, el("\"Show side chains\" checkbox")));
      await session.step(112, "Then \"Show side chains\" checkbox should be checked", () => shouldBe(page, el("\"Show side chains\" checkbox"), "checked"));
      await session.step(113, "And \"showBindingSite\" property of Biostructure viewer should be \"true\"", () => propertyShouldBe(page, "showBindingSite", el("Biostructure viewer"), "true"));
      await session.step(114, "And \"Binding site\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Binding site\" overlay button in Biostructure viewer"), "selected"));
      await session.step(115, "And the \"binding site shown\" reading of Biostructure viewer should be \"true\"", () => readingReads(page, "binding site shown", el("Biostructure viewer"), "true"));
      await session.step(116, "And the \"binding site atoms\" reading of Biostructure viewer should be 335", () => readingIs(page, "binding site atoms", el("Biostructure viewer"), 335));
      await session.step(117, "Given \"Binding Site\" category in context panel is expanded", () => isExpanded(page, el("\"Binding Site\" category in context panel")));
      await session.step(118, "Then \"Show Binding Site\" property in context panel should be checked", () => shouldBe(page, el("\"Show Binding Site\" property in context panel"), "checked"));
      await session.step(119, "When user enters \"8\" in \"Binding Site Radius\" property in context panel", () => enterInto(page, "8", el("\"Binding Site Radius\" property in context panel")));
      await session.step(120, "Then \"bindingSiteRadius\" property of Biostructure viewer should be \"8\"", () => propertyShouldBe(page, "bindingSiteRadius", el("Biostructure viewer"), "8"));
      await session.step(121, "And the \"binding site atoms\" reading of Biostructure viewer should be 559", () => readingIs(page, "binding site atoms", el("Biostructure viewer"), 559));
      await session.step(122, "When user clicks on \"Binding site\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Binding site\" overlay button in Biostructure viewer")));
      await session.step(123, "Then \"Show side chains\" checkbox should be visible", () => shouldBe(page, el("\"Show side chains\" checkbox"), "visible"));
      await session.step(124, "And \"8.0 Å\" text should be visible", () => shouldBe(page, el("\"8.0 Å\" text"), "visible"));
      await session.step(125, "When user unchecks \"Show Binding Site\" property in context panel", () => uncheck(page, el("\"Show Binding Site\" property in context panel")));
      await session.step(126, "Then \"showBindingSite\" property of Biostructure viewer should be \"false\"", () => propertyShouldBe(page, "showBindingSite", el("Biostructure viewer"), "false"));
      await session.step(127, "And \"Binding site\" overlay button in Biostructure viewer should not be selected", () => shouldNotBe(page, el("\"Binding site\" overlay button in Biostructure viewer"), "selected"));
      await session.step(128, "And the \"binding site shown\" reading of Biostructure viewer should be \"false\"", () => readingReads(page, "binding site shown", el("Biostructure viewer"), "false"));
      await session.step(129, "When user clicks on \"Binding site\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Binding site\" overlay button in Biostructure viewer")));
      await session.step(130, "Then \"Show side chains\" checkbox should be visible", () => shouldBe(page, el("\"Show side chains\" checkbox"), "visible"));
      await session.step(131, "And \"Show side chains\" checkbox should be unchecked", () => shouldBe(page, el("\"Show side chains\" checkbox"), "unchecked"));
      await session.step(132, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(133, "Then \"Show side chains\" checkbox should be hidden", () => shouldBe(page, el("\"Show side chains\" checkbox"), "hidden"));
      await session.step(134, "When user enters \"5\" in \"Binding Site Radius\" property in context panel", () => enterInto(page, "5", el("\"Binding Site Radius\" property in context panel")));
      await session.step(135, "Then \"bindingSiteRadius\" property of Biostructure viewer should be \"5\"", () => propertyShouldBe(page, "bindingSiteRadius", el("Biostructure viewer"), "5"));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toggle Expanded Viewport expands the Mol* layout and Escape brings it back", async () => {
      await session.step(139, "Then \"Assembly\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "absent"));
      await session.step(140, "And the \"layout expanded\" reading of Biostructure viewer should be \"false\"", () => readingReads(page, "layout expanded", el("Biostructure viewer"), "false"));
      await session.step(141, "When user clicks on \"Toggle Expanded Viewport\" overlay button in Biostructure viewer", () => clickOn(page, el("\"Toggle Expanded Viewport\" overlay button in Biostructure viewer")));
      await session.step(142, "Then the \"layout expanded\" reading of Biostructure viewer should be \"true\"", () => readingReads(page, "layout expanded", el("Biostructure viewer"), "true"));
      await session.step(143, "And \"Toggle Expanded Viewport\" overlay button in Biostructure viewer should be selected", () => shouldBe(page, el("\"Toggle Expanded Viewport\" overlay button in Biostructure viewer"), "selected"));
      await session.step(144, "And \"Assembly\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "visible"));
      await session.step(145, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(146, "Then the \"layout expanded\" reading of Biostructure viewer should be \"false\"", () => readingReads(page, "layout expanded", el("Biostructure viewer"), "false"));
      await session.step(147, "And \"Assembly\" button in Biostructure viewer should be absent", () => shouldBe(page, el("\"Assembly\" button in Biostructure viewer"), "absent"));
      await session.step(148, "And the \"pdb_data\" view should be current", () => viewIsCurrent(page, "pdb_data"));
      await session.step(149, "And Biostructure viewer should be visible", () => shouldBe(page, el("Biostructure viewer"), "visible"));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
