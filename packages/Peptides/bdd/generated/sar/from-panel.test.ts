/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/from-panel.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, expand, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, painted, readingIs, readingNotAsRemembered, readingReads, rememberReading, repainted, takeSnapshot, viewerAdded, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Launch and configure SAR from the Peptides pane", () => {
  const session = feature(test, "features/sar/from-panel.feature", import.meta.url);
  test("Launch and configure SAR from the Peptides pane", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(15, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "When user clicks on the \"header AlignedSequence\" area of grid", () => clickArea(page, "header AlignedSequence", el("grid")));
    await session.step(18, "Then the context panel should show \"AlignedSequence\"", () => contextPanelShows(page, "AlignedSequence"));
    await session.step(19, "When user expands Peptides pane in context panel", () => expand(page, el("Peptides pane in context panel")));
    await session.step(20, "Then \"Launch SAR\" button in Peptides pane should be visible", () => shouldBe(page, el("\"Launch SAR\" button in Peptides pane"), "visible"));
    await session.step(21, "When user clicks on \"Adjust clustering parameters\" icon in Peptides pane", () => clickOn(page, el("\"Adjust clustering parameters\" icon in Peptides pane")));
    await session.step(22, "And user enters \"93\" into \"Similarity Threshold\" input in Peptides pane", () => enterInto(page, "93", el("\"Similarity Threshold\" input in Peptides pane")));
    await session.step(23, "Then no errors should have been logged", () => noErrors(page));
    await session.step(24, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Launch SAR attaches the analysis and clustering results", async () => {
      await session.step(27, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(28, "When user clicks on \"Launch SAR\" button in Peptides pane", () => clickOn(page, el("\"Launch SAR\" button in Peptides pane")));
      await session.step(29, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(30, "And Sequence Variability Map viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Variability Map"));
      await session.step(31, "And Most Potent Residues viewer should be added to the open tableview", () => viewerAdded(page, "Most Potent Residues"));
      await session.step(32, "And MCL viewer should be added to the open tableview", () => viewerAdded(page, "MCL"));
      await session.step(33, "And Logo Summary Table viewer should be added to the open tableview", () => viewerAdded(page, "Logo Summary Table"));
      await session.step(34, "And the table should have a column \"Cluster (MCL)\"", () => hasColumn(page, "Cluster (MCL)"));
      await session.step(35, "And the \"clusters column\" reading of Logo Summary Table viewer should be \"Cluster (MCL)\"", () => readingReads(page, "clusters column", el("Logo Summary Table viewer"), "Cluster (MCL)"));
      await session.step(36, "And the \"members total\" reading of Logo Summary Table viewer should be 647", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 647));
      await session.step(37, "And the \"completed threshold\" reading of MCL viewer should be 93", () => readingIs(page, "completed threshold", el("MCL viewer"), 93));
      await session.step(38, "And Sequence Variability Map viewer should be painted", () => painted(page, el("Sequence Variability Map viewer")));
      await session.step(39, "And Most Potent Residues viewer should be painted", () => painted(page, el("Most Potent Residues viewer")));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The settings wrench exposes the launch MCL configuration", async () => {
      await session.step(44, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(45, "Then \"Peptides settings\" dialog should be visible", () => shouldBe(page, el("\"Peptides settings\" dialog"), "visible"));
      await session.step(46, "When user expands MCL pane in \"Peptides settings\" dialog", () => expand(page, el("MCL pane in \"Peptides settings\" dialog")));
      await session.step(47, "Then \"Similarity Threshold\" input in \"Peptides settings\" dialog should have value \"93\"", () => shouldHaveValue(page, el("\"Similarity Threshold\" input in \"Peptides settings\" dialog"), "93"));
      await session.step(48, "And \"Inflation Factor\" input in \"Peptides settings\" dialog should have value \"1.4\"", () => shouldHaveValue(page, el("\"Inflation Factor\" input in \"Peptides settings\" dialog"), "1.4"));
      await session.step(49, "When user clicks on CANCEL button in \"Peptides settings\" dialog", () => clickOn(page, el("CANCEL button in \"Peptides settings\" dialog")));
      await session.step(50, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(51, "Then no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A changed similarity threshold completes clustering and preserves the sequence viewers", async () => {
      await session.step(55, "When user remembers the \"completed computations\" reading of MCL viewer", () => rememberReading(page, "completed computations", el("MCL viewer")));
      await session.step(56, "And user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(57, "And user expands MCL pane in \"Peptides settings\" dialog", () => expand(page, el("MCL pane in \"Peptides settings\" dialog")));
      await session.step(58, "And user enters \"90\" into \"Similarity Threshold\" input in \"Peptides settings\" dialog", () => enterInto(page, "90", el("\"Similarity Threshold\" input in \"Peptides settings\" dialog")));
      await session.step(59, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(60, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(61, "Then \"Peptides settings\" dialog should be hidden", () => shouldBe(page, el("\"Peptides settings\" dialog"), "hidden"));
      await session.step(62, "And the SAR analysis should be ready", () => sarReady(page));
      await session.step(63, "And the \"completed computations\" reading of MCL viewer should not be as remembered", () => readingNotAsRemembered(page, "completed computations", el("MCL viewer")));
      await session.step(64, "And the SAR setting \"mclSettings.threshold\" should be \"90\"", () => sarSetting(page, "mclSettings.threshold", "90"));
      await session.step(65, "And the \"completed threshold\" reading of MCL viewer should be 90", () => readingIs(page, "completed threshold", el("MCL viewer"), 90));
      await session.step(66, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(67, "And scatter plot viewer in MCL viewer should be painted", () => painted(page, el("scatter plot viewer in MCL viewer")));
      await session.step(68, "And Sequence Variability Map viewer should be visible", () => shouldBe(page, el("Sequence Variability Map viewer"), "visible"));
      await session.step(69, "And Most Potent Residues viewer should be visible", () => shouldBe(page, el("Most Potent Residues viewer"), "visible"));
      await session.step(70, "And the \"members total\" reading of Logo Summary Table viewer should be 647", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 647));
      await session.step(71, "When user remembers the \"completed computations\" reading of MCL viewer", () => rememberReading(page, "completed computations", el("MCL viewer")));
      await session.step(72, "And user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(73, "And user expands MCL pane in \"Peptides settings\" dialog", () => expand(page, el("MCL pane in \"Peptides settings\" dialog")));
      await session.step(74, "And user enters \"93\" into \"Similarity Threshold\" input in \"Peptides settings\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Peptides settings\" dialog")));
      await session.step(75, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(76, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(77, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(78, "And the \"completed computations\" reading of MCL viewer should not be as remembered", () => readingNotAsRemembered(page, "completed computations", el("MCL viewer")));
      await session.step(79, "And the SAR setting \"mclSettings.threshold\" should be \"93\"", () => sarSetting(page, "mclSettings.threshold", "93"));
      await session.step(80, "And the \"completed threshold\" reading of MCL viewer should be 93", () => readingIs(page, "completed threshold", el("MCL viewer"), 93));
      await session.step(81, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The invariant map selects peptides and opens their position distribution", async () => {
      await session.step(86, "When user takes a snapshot of Sequence Variability Map viewer", () => takeSnapshot(page, el("Sequence Variability Map viewer")));
      await session.step(87, "And user checks \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => check(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(88, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Invariant Map\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Invariant Map"));
      await session.step(89, "And Sequence Variability Map viewer should have repainted", () => repainted(page, el("Sequence Variability Map viewer")));
      await session.step(90, "And \"Mutation Cliffs\" checkbox in Sequence Variability Map viewer should be unchecked", () => shouldBe(page, el("\"Mutation Cliffs\" checkbox in Sequence Variability Map viewer"), "unchecked"));
      await session.step(91, "When user clicks on the \"cell A at 2\" area of Sequence Variability Map viewer", () => clickArea(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(92, "Then only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
      await session.step(93, "And 299 rows should be selected", () => selectedRowCount(page, 299));
      await session.step(94, "And the \"selected monomer-positions\" reading of Sequence Variability Map viewer should be \"2:A\"", () => readingReads(page, "selected monomer-positions", el("Sequence Variability Map viewer"), "2:A"));
      await session.step(95, "And \"Mutation Cliffs pairs\" pane in context panel should be present", () => shouldBe(page, el("\"Mutation Cliffs pairs\" pane in context panel"), "present"));
      await session.step(96, "And Distribution pane in context panel should be visible", () => shouldBe(page, el("Distribution pane in context panel"), "visible"));
      await session.step(97, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(98, "Then Distribution pane in context panel should contain text \"Mean difference\"", () => shouldContainText(page, el("Distribution pane in context panel"), "Mean difference"));
      await session.step(99, "And Distribution pane in context panel should not contain text \"No distribution\"", () => shouldNotContainText(page, el("Distribution pane in context panel"), "No distribution"));
      await session.step(100, "And Positions heading in Distribution pane should be hidden", () => shouldBe(page, el("Positions heading in Distribution pane"), "hidden"));
      await session.step(101, "When user checks Positions checkbox in Distribution pane", () => check(page, el("Positions checkbox in Distribution pane")));
      await session.step(102, "Then Positions heading in Distribution pane should be visible", () => shouldBe(page, el("Positions heading in Distribution pane"), "visible"));
      await session.step(103, "And second Histogram viewer in Distribution pane should be painted", () => painted(page, el("second Histogram viewer in Distribution pane")));
      await session.step(104, "When user unchecks Positions checkbox in Distribution pane", () => uncheck(page, el("Positions checkbox in Distribution pane")));
      await session.step(105, "Then Positions heading in Distribution pane should be hidden", () => shouldBe(page, el("Positions heading in Distribution pane"), "hidden"));
      await session.step(106, "When user clears the row selection", () => clearSelection(page));
      await session.step(107, "And user takes a snapshot of Sequence Variability Map viewer", () => takeSnapshot(page, el("Sequence Variability Map viewer")));
      await session.step(108, "And user checks \"Mutation Cliffs\" checkbox in Sequence Variability Map viewer", () => check(page, el("\"Mutation Cliffs\" checkbox in Sequence Variability Map viewer")));
      await session.step(109, "Then no rows should be selected", () => noneSelected(page));
      await session.step(110, "And the \"mode\" reading of Sequence Variability Map viewer should be \"Mutation Cliffs\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Mutation Cliffs"));
      await session.step(111, "And Sequence Variability Map viewer should have repainted", () => repainted(page, el("Sequence Variability Map viewer")));
      await session.step(112, "And \"Invariant Map\" checkbox in Sequence Variability Map viewer should be unchecked", () => shouldBe(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer"), "unchecked"));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
      await session.step(114, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
