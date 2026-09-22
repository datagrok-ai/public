/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/tooltips.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clusterStatistics, clusterSummaries} from '../../bindings/clusters.js';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, expand, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {contextPanelOpen, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noBalloons, noErrors, pointerAway, readingIs, readingReads, repainted, takeSnapshot, viewerAdded, viewerCount, wheelOverAreaTimesHolding} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Inspect peptide statistics in tooltips", () => {
  const session = feature(test, "features/sar/tooltips.feature", import.meta.url);
  test("Inspect peptide statistics in tooltips", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(9, "And user opens peptides dataset keeping the first 100 rows", () => openDatasetRows(page, ds("peptides"), 100));
    await session.step(10, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(11, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(12, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(13, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(14, "And user enters \"94\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "94", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(15, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(16, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(17, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(18, "And Sequence Variability Map viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Variability Map"));
    await session.step(19, "And Most Potent Residues viewer should be added to the open tableview", () => viewerAdded(page, "Most Potent Residues"));
    await session.step(20, "Then no errors should have been logged", () => noErrors(page));
    await session.step(21, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Invariant-map tooltips show the hovered cell's population", async () => {
      await session.step(24, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(25, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Invariant Map\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Invariant Map"));
      await session.step(26, "When user hovers over the \"cell A at 2\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(27, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(28, "And Count table row in tooltip should contain text \"14 (14.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "14 (14.000%)"));
      await session.step(29, "And \"Mean difference\" table row in tooltip should be visible", () => shouldBe(page, el("\"Mean difference\" table row in tooltip"), "visible"));
      await session.step(30, "And the \"count of cell A at 2\" reading of Sequence Variability Map viewer should be 14", () => readingIs(page, "count of cell A at 2", el("Sequence Variability Map viewer"), 14));
      await session.step(31, "And the \"highlighted rows\" reading of grid should be 14", () => readingIs(page, "highlighted rows", el("grid"), 14));
      await session.step(32, "When user hovers over the \"cell N at 4\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell N at 4", el("Sequence Variability Map viewer")));
      await session.step(33, "Then Count table row in tooltip should contain text \"74 (74.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "74 (74.000%)"));
      await session.step(34, "And Count table row in tooltip should not contain text \"14 (14.000%)\"", () => shouldNotContainText(page, el("Count table row in tooltip"), "14 (14.000%)"));
      await session.step(35, "And the \"highlighted rows\" reading of grid should be 74", () => readingIs(page, "highlighted rows", el("grid"), 74));
      await session.step(36, "Then no errors should have been logged", () => noErrors(page));
      await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Mutation-cliff tooltips report the number of substitution pairs", async () => {
      await session.step(40, "When user takes a snapshot of Sequence Variability Map viewer", () => takeSnapshot(page, el("Sequence Variability Map viewer")));
      await session.step(41, "When user clicks on \"Mutation Cliffs\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Mutation Cliffs\" checkbox in Sequence Variability Map viewer")));
      await session.step(42, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Mutation Cliffs\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Mutation Cliffs"));
      await session.step(43, "And Sequence Variability Map viewer should have repainted", () => repainted(page, el("Sequence Variability Map viewer")));
      await session.step(44, "And the \"cliffs of cell A at 2\" reading of Sequence Variability Map viewer should be 7", () => readingIs(page, "cliffs of cell A at 2", el("Sequence Variability Map viewer"), 7));
      await session.step(45, "When user hovers over the \"cell A at 2\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(46, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(47, "And \"Pairs count\" table row in tooltip should have text \"Pairs count7\"", () => shouldHaveText(page, el("\"Pairs count\" table row in tooltip"), "Pairs count7"));
      await session.step(48, "And \"MP in cliffs count\" table row in tooltip should be visible", () => shouldBe(page, el("\"MP in cliffs count\" table row in tooltip"), "visible"));
      await session.step(49, "Then no errors should have been logged", () => noErrors(page));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A WebLogo selection leaves the invariant-map tooltip usable", async () => {
      await session.step(53, "When user moves the pointer away from Sequence Variability Map viewer", () => pointerAway(page, el("Sequence Variability Map viewer")));
      await session.step(54, "And user clicks on the \"A at 2\" area of grid", () => clickArea(page, "A at 2", el("grid")));
      await session.step(55, "Then 14 rows should be selected", () => selectedRowCount(page, 14));
      await session.step(56, "And only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
      await session.step(57, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(58, "And user hovers over the \"cell N at 4\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell N at 4", el("Sequence Variability Map viewer")));
      await session.step(59, "Then Count table row in tooltip should contain text \"74 (74.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "74 (74.000%)"));
      await session.step(60, "And the \"highlighted rows\" reading of grid should be 74", () => readingIs(page, "highlighted rows", el("grid"), 74));
      await session.step(61, "When user moves the pointer away from Sequence Variability Map viewer", () => pointerAway(page, el("Sequence Variability Map viewer")));
      await session.step(62, "Then Count table row in tooltip should be hidden", () => shouldBe(page, el("Count table row in tooltip"), "hidden"));
      await session.step(63, "And the \"highlighted rows\" reading of grid should be 0", () => readingIs(page, "highlighted rows", el("grid"), 0));
      await session.step(64, "Then no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Adding and removing the active-peptide viewer preserves tooltips", async () => {
      await session.step(68, "When user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(69, "Then \"Peptides settings\" dialog should be visible", () => shouldBe(page, el("\"Peptides settings\" dialog"), "visible"));
      await session.step(70, "When user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(71, "And user checks \"Active peptide selection\" checkbox in \"Peptides settings\" dialog", () => check(page, el("\"Active peptide selection\" checkbox in \"Peptides settings\" dialog")));
      await session.step(72, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(73, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(74, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(75, "And Active peptide selection viewer should be added to the open tableview", () => viewerAdded(page, "Active peptide selection"));
      await session.step(76, "When user hovers over the \"cell A at 2\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(77, "Then Count table row in tooltip should contain text \"14 (14.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "14 (14.000%)"));
      await session.step(78, "When user moves the pointer away from Sequence Variability Map viewer", () => pointerAway(page, el("Sequence Variability Map viewer")));
      await session.step(79, "And user clicks on \"Peptides analysis settings\" icon", () => clickOn(page, el("\"Peptides analysis settings\" icon")));
      await session.step(80, "And user expands Viewers pane in \"Peptides settings\" dialog", () => expand(page, el("Viewers pane in \"Peptides settings\" dialog")));
      await session.step(81, "And user unchecks \"Active peptide selection\" checkbox in \"Peptides settings\" dialog", () => uncheck(page, el("\"Active peptide selection\" checkbox in \"Peptides settings\" dialog")));
      await session.step(82, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(83, "When user clicks on OK button in \"Peptides settings\" dialog", () => clickOn(page, el("OK button in \"Peptides settings\" dialog")));
      await session.step(84, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(85, "And the open tableview should have 0 Active peptide selection viewers", () => viewerCount(page, 0, "Active peptide selection"));
      await session.step(86, "When user hovers over the \"cell N at 4\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell N at 4", el("Sequence Variability Map viewer")));
      await session.step(87, "Then Count table row in tooltip should contain text \"74 (74.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "74 (74.000%)"));
      await session.step(88, "Then no errors should have been logged", () => noErrors(page));
      await session.step(89, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("WebLogo header tooltips replace their statistics and clear the highlight on leaving", async () => {
      await session.step(94, "When user moves the pointer away from Sequence Variability Map viewer", () => pointerAway(page, el("Sequence Variability Map viewer")));
      await session.step(95, "And user scrolls the mouse wheel up 5 times over the \"row header 1\" area of grid holding Shift", () => wheelOverAreaTimesHolding(page, "up", 5, "row header 1", el("grid"), "Shift"));
      await session.step(96, "And user hovers over the \"A at 2\" area of grid", () => hoverArea(page, "A at 2", el("grid")));
      await session.step(97, "Then Count table row in tooltip should contain text \"14 (14.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "14 (14.000%)"));
      await session.step(98, "And the \"highlighted rows\" reading of grid should be 14", () => readingIs(page, "highlighted rows", el("grid"), 14));
      await session.step(99, "When user hovers over the \"N at 4\" area of grid", () => hoverArea(page, "N at 4", el("grid")));
      await session.step(100, "Then Count table row in tooltip should contain text \"74 (74.000%)\"", () => shouldContainText(page, el("Count table row in tooltip"), "74 (74.000%)"));
      await session.step(101, "And the \"highlighted rows\" reading of grid should be 74", () => readingIs(page, "highlighted rows", el("grid"), 74));
      await session.step(102, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(103, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(104, "And the \"highlighted rows\" reading of grid should be 0", () => readingIs(page, "highlighted rows", el("grid"), 0));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
      await session.step(106, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Cluster tooltips and selection use the computed cluster membership", async () => {
      await session.step(109, "When user clears the row selection", () => clearSelection(page));
      await session.step(110, "Then no rows should be selected", () => noneSelected(page));
      await session.step(111, "And the \"members total\" reading of Logo Summary Table viewer should be 100", () => readingIs(page, "members total", el("Logo Summary Table viewer"), 100));
      await session.step(112, "And the cluster summaries should match the source rows of table \"peptides\"", () => clusterSummaries(page, "peptides"));
      await session.step(113, "When user hovers over the \"cluster 1\" area of Logo Summary Table viewer", () => hoverArea(page, "cluster 1", el("Logo Summary Table viewer")));
      await session.step(114, "Then the cluster statistics in tooltip should match cluster \"1\" of table \"peptides\"", () => clusterStatistics(page, el("tooltip"), "1", "peptides"));
      await session.step(115, "When user hovers over the \"cluster 2\" area of Logo Summary Table viewer", () => hoverArea(page, "cluster 2", el("Logo Summary Table viewer")));
      await session.step(116, "Then the cluster statistics in tooltip should match cluster \"2\" of table \"peptides\"", () => clusterStatistics(page, el("tooltip"), "2", "peptides"));
      await session.step(117, "When user clicks on the \"cluster 1\" area of Logo Summary Table viewer", () => clickArea(page, "cluster 1", el("Logo Summary Table viewer")));
      await session.step(118, "Then only rows where \"Cluster (MCL)\" is \"1\" should be selected", () => onlyOfSelected(page, "Cluster (MCL)", "1"));
      await session.step(119, "And the \"selected clusters\" reading of Logo Summary Table viewer should be \"1\"", () => readingReads(page, "selected clusters", el("Logo Summary Table viewer"), "1"));
      await session.step(120, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(121, "Then the cluster statistics in Distribution pane in context panel should match cluster \"1\" of table \"peptides\"", () => clusterStatistics(page, el("Distribution pane in context panel"), "1", "peptides"));
      await session.step(122, "When user clicks on the \"cluster 2\" area of Logo Summary Table viewer", () => clickArea(page, "cluster 2", el("Logo Summary Table viewer")));
      await session.step(123, "Then only rows where \"Cluster (MCL)\" is \"2\" should be selected", () => onlyOfSelected(page, "Cluster (MCL)", "2"));
      await session.step(124, "And the \"selected clusters\" reading of Logo Summary Table viewer should be \"2\"", () => readingReads(page, "selected clusters", el("Logo Summary Table viewer"), "2"));
      await session.step(125, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(126, "Then the cluster statistics in Distribution pane in context panel should match cluster \"2\" of table \"peptides\"", () => clusterStatistics(page, el("Distribution pane in context panel"), "2", "peptides"));
      await session.step(127, "When user clears the row selection", () => clearSelection(page));
      await session.step(128, "Then no rows should be selected", () => noneSelected(page));
      await session.step(129, "And the \"selected clusters\" reading of Logo Summary Table viewer should be \"\"", () => readingReads(page, "selected clusters", el("Logo Summary Table viewer"), ""));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
      await session.step(131, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
