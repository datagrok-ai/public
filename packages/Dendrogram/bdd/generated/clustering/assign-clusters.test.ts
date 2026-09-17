/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/clustering/assign-clusters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [dendrogram.cp.assign-clusters-column-creation]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldContainText, shouldHaveValue, shouldHaveValueBetween, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnType, currentRowIs, hasColumn, mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnsMatching, newestMatchingDistinct, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterBetween, filterPassesFewer, resetFilter, selectFirstRows, selectNoRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hasArea, hasNoArea, hoverArea, noBalloons, noErrors, pickFromAreaContextMenu, readingHigher, readingIs, readingLower, readingReads, readingsEqual, warningBalloonText, wheelOverArea, wheelOverAreaHolding, wheelOverAreaTimesHolding} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Assigning clusters from the tree next to the grid", () => {
  const session = feature(test, "features/clustering/assign-clusters.feature", import.meta.url);
  test("Assigning clusters from the tree next to the grid", {tag: ["@journey", "@realizes:dendrogram.cp.assign-clusters-column-creation", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens mol1K dataset", () => openDataset(page, ds("mol1K")));
    await session.step(23, "And user watches the task bar", () => watchTaskBar(page));
    await session.step(24, "When user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
    await session.step(25, "And user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
    await session.step(26, "Then the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
    await run.scenario("The grid draws structures and the clustering shows its progress", async () => {
      await session.step(29, "Then the \"cell type of molecule\" reading of grid should be \"Molecule\"", () => readingReads(page, "cell type of molecule", el("grid"), "Molecule"));
      await session.step(30, "And the task bar should have shown \"Creating dendrogram\"", () => taskBarShown(page, "Creating dendrogram"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The tree follows the grid's current row, selection and filter", async () => {
      await session.step(34, "When user clicks on the \"cell 3 of prID\" area of grid", () => clickArea(page, "cell 3 of prID", el("grid")));
      await session.step(35, "Then row 3 should be current", () => currentRowIs(page, 3));
      await session.step(36, "And the \"tree current node\" reading of grid should be \"2\"", () => readingReads(page, "tree current node", el("grid"), "2"));
      await session.step(37, "When user selects the first 10 rows", () => selectFirstRows(page, 10));
      await session.step(38, "Then the \"tree selected leaves\" reading of grid should be 10", () => readingIs(page, "tree selected leaves", el("grid"), 10));
      await session.step(39, "When user filters rows where \"pIC50_HIV_Integrase\" is between 6 and 7", () => filterBetween(page, "pIC50_HIV_Integrase", 6, 7));
      await session.step(40, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(41, "And the \"tree leaves\" reading of grid should be lower than before", () => readingLower(page, "tree leaves", el("grid")));
      await session.step(42, "And the \"tree leaves\" and \"rows shown\" readings of grid should be the same", () => readingsEqual(page, "tree leaves", "rows shown", el("grid")));
      await session.step(43, "When user resets the filter", () => resetFilter(page));
      await session.step(44, "And user selects no rows", () => selectNoRows(page));
      await session.step(45, "Then the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(46, "And the \"tree selected leaves\" reading of grid should be 0", () => readingIs(page, "tree selected leaves", el("grid"), 0));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid follows the tree: a leaf click, a hover, an inner node click", async () => {
      await session.step(50, "When user clicks on the \"cell 266 of prID\" area of grid", () => clickArea(page, "cell 266 of prID", el("grid")));
      await session.step(51, "Then the \"tree current node\" reading of grid should be \"265\"", () => readingReads(page, "tree current node", el("grid"), "265"));
      await session.step(52, "When user clicks on the \"leaf 14\" area of grid", () => clickArea(page, "leaf 14", el("grid")));
      await session.step(53, "Then row 15 should be current", () => currentRowIs(page, 15));
      await session.step(54, "And the \"tree current node\" reading of grid should be \"14\"", () => readingReads(page, "tree current node", el("grid"), "14"));
      await session.step(55, "When user hovers over the \"leaf 104\" area of grid", () => hoverArea(page, "leaf 104", el("grid")));
      await session.step(56, "Then the \"tree mouse over node\" reading of grid should be \"104\"", () => readingReads(page, "tree mouse over node", el("grid"), "104"));
      await session.step(57, "And row 105 should be under the mouse", () => mouseOverRowIs(page, 105));
      await session.step(58, "When user clicks on the \"node 265-393\" area of grid", () => clickArea(page, "node 265-393", el("grid")));
      await session.step(59, "Then 41 rows should be selected", () => selectedRowCount(page, 41));
      await session.step(60, "And the \"tree selected leaves\" reading of grid should be 41", () => readingIs(page, "tree selected leaves", el("grid"), 41));
      await session.step(61, "When user selects no rows", () => selectNoRows(page));
      await session.step(62, "Then the \"tree selected leaves\" reading of grid should be 0", () => readingIs(page, "tree selected leaves", el("grid"), 0));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The tree's context menu opens Assign Clusters at half the tree's height, with its cut line there", async () => {
      await session.step(66, "When user picks \"Assign Clusters\" from the context menu of the \"tree\" area of grid", () => pickFromAreaContextMenu(page, "Assign Clusters", "tree", el("grid")));
      await session.step(67, "Then \"Assign Clusters\" dialog should be visible", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "visible"));
      await session.step(68, "And Threshold input in \"Assign Clusters\" dialog should have a value between 319.5 and 319.6", () => shouldHaveValueBetween(page, el("Threshold input in \"Assign Clusters\" dialog"), 319.5, 319.6));
      await session.step(69, "And Clusters input in \"Assign Clusters\" dialog should have value \"6\"", () => shouldHaveValue(page, el("Clusters input in \"Assign Clusters\" dialog"), "6"));
      await session.step(70, "And \"Medoid columns\" input in \"Assign Clusters\" dialog should be checked", () => shouldBe(page, el("\"Medoid columns\" input in \"Assign Clusters\" dialog"), "checked"));
      await session.step(71, "And the \"cut threshold\" reading of grid should be 319.54", () => readingIs(page, "cut threshold", el("grid"), 319.54));
      await session.step(72, "And grid should have a \"cut line\" area", () => hasArea(page, el("grid"), "cut line"));
      await session.step(73, "When user clicks on CANCEL button in \"Assign Clusters\" dialog", () => clickOn(page, el("CANCEL button in \"Assign Clusters\" dialog")));
      await session.step(74, "Then \"Assign Clusters\" dialog should be hidden", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "hidden"));
      await session.step(75, "And grid should not have a \"cut line\" area", () => hasNoArea(page, el("grid"), "cut line"));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The magic wand opens the same dialog", async () => {
      await session.step(79, "When user hovers over \"Assign Clusters\" icon", () => hoverOver(page, el("\"Assign Clusters\" icon")));
      await session.step(80, "Then tooltip should contain text \"Assign Clusters\"", () => shouldContainText(page, el("tooltip"), "Assign Clusters"));
      await session.step(81, "When user clicks on \"Assign Clusters\" icon", () => clickOn(page, el("\"Assign Clusters\" icon")));
      await session.step(82, "Then \"Assign Clusters\" dialog should be visible", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "visible"));
      await session.step(83, "And Threshold input in \"Assign Clusters\" dialog should have a value between 319.5 and 319.6", () => shouldHaveValueBetween(page, el("Threshold input in \"Assign Clusters\" dialog"), 319.5, 319.6));
      await session.step(84, "When user clicks on CANCEL button in \"Assign Clusters\" dialog", () => clickOn(page, el("CANCEL button in \"Assign Clusters\" dialog")));
      await session.step(85, "Then \"Assign Clusters\" dialog should be hidden", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "hidden"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Threshold sets Clusters and moves the cut line, Clusters sets Threshold, and Assign adds the columns of that cut", async () => {
      await session.step(89, "When user clicks on \"Assign Clusters\" icon", () => clickOn(page, el("\"Assign Clusters\" icon")));
      await session.step(90, "And user enters \"200\" into Threshold input in \"Assign Clusters\" dialog", () => enterInto(page, "200", el("Threshold input in \"Assign Clusters\" dialog")));
      await session.step(91, "Then Clusters input in \"Assign Clusters\" dialog should have value \"4\"", () => shouldHaveValue(page, el("Clusters input in \"Assign Clusters\" dialog"), "4"));
      await session.step(92, "And the \"cut threshold\" reading of grid should be 200", () => readingIs(page, "cut threshold", el("grid"), 200));
      await session.step(93, "When user enters \"5\" into Clusters input in \"Assign Clusters\" dialog", () => enterInto(page, "5", el("Clusters input in \"Assign Clusters\" dialog")));
      await session.step(94, "Then Threshold input in \"Assign Clusters\" dialog should have a value between 239.65 and 239.67", () => shouldHaveValueBetween(page, el("Threshold input in \"Assign Clusters\" dialog"), 239.65, 239.67));
      await session.step(95, "And the \"cut threshold\" reading of grid should be 239.66", () => readingIs(page, "cut threshold", el("grid"), 239.66));
      await session.step(96, "When user clicks on Assign button in \"Assign Clusters\" dialog", () => clickOn(page, el("Assign button in \"Assign Clusters\" dialog")));
      await session.step(97, "Then the \"Assign Clusters\" dialog should close", () => dialogCloses(page, "Assign Clusters"));
      await session.step(98, "And 3 new columns matching \"\\(239\\.66\\)$\" should have been added", () => newColumnsMatching(page, 3, "\\(239\\.66\\)$"));
      await session.step(99, "And the table should have a column \"Cluster (239.66)\"", () => hasColumn(page, "Cluster (239.66)"));
      await session.step(100, "And the table should have a column \"Medoid Rank (239.66)\"", () => hasColumn(page, "Medoid Rank (239.66)"));
      await session.step(101, "And the table should have a column \"Avg Distance to Cluster (239.66)\"", () => hasColumn(page, "Avg Distance to Cluster (239.66)"));
      await session.step(102, "And \"Cluster (239.66)\" column should have type \"string\"", () => columnType(page, "Cluster (239.66)", "string"));
      await session.step(103, "And \"Cluster (239.66)\" column should have no missing values", () => columnComplete(page, "Cluster (239.66)"));
      await session.step(104, "And the newest column matching \"^Cluster \\(\" should have 5 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(", 5));
      await session.step(105, "And \"Medoid Rank (239.66)\" column should have no missing values", () => columnComplete(page, "Medoid Rank (239.66)"));
      await session.step(106, "And \"Avg Distance to Cluster (239.66)\" column should have no missing values", () => columnComplete(page, "Avg Distance to Cluster (239.66)"));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second Assign adds its own columns and keeps the first", async () => {
      await session.step(111, "When user clicks on \"Assign Clusters\" icon", () => clickOn(page, el("\"Assign Clusters\" icon")));
      await session.step(112, "And user enters \"3\" into Clusters input in \"Assign Clusters\" dialog", () => enterInto(page, "3", el("Clusters input in \"Assign Clusters\" dialog")));
      await session.step(113, "And user clicks on Assign button in \"Assign Clusters\" dialog", () => clickOn(page, el("Assign button in \"Assign Clusters\" dialog")));
      await session.step(114, "Then the \"Assign Clusters\" dialog should close", () => dialogCloses(page, "Assign Clusters"));
      await session.step(115, "And 3 new columns matching \"\\(159\\.77\\)$\" should have been added", () => newColumnsMatching(page, 3, "\\(159\\.77\\)$"));
      await session.step(116, "And 2 new columns matching \"^Cluster \\(\" should have been added", () => newColumnsMatching(page, 2, "^Cluster \\("));
      await session.step(117, "And the table should have a column \"Cluster (239.66)\"", () => hasColumn(page, "Cluster (239.66)"));
      await session.step(118, "And the newest column matching \"^Cluster \\(\" should have 3 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(", 3));
      await session.step(119, "And the newest column matching \"^Cluster \\(\" should have no missing values", () => newestMatchingFilled(page, "^Cluster \\("));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wheel scrolls the tree with the grid, and Control with the wheel zooms it between 1 and 100", async () => {
      await session.step(123, "When user scrolls the mouse wheel down over the \"tree\" area of grid", () => wheelOverArea(page, "down", "tree", el("grid")));
      await session.step(124, "Then the \"tree top row\" reading of grid should be higher than before", () => readingHigher(page, "tree top row", el("grid")));
      await session.step(125, "When user scrolls the mouse wheel up over the \"tree\" area of grid", () => wheelOverArea(page, "up", "tree", el("grid")));
      await session.step(126, "Then the \"tree top row\" reading of grid should be 0", () => readingIs(page, "tree top row", el("grid"), 0));
      await session.step(127, "When user scrolls the mouse wheel up over the \"tree\" area of grid holding Control", () => wheelOverAreaHolding(page, "up", "tree", el("grid"), "Control"));
      await session.step(128, "Then the \"tree zoom\" reading of grid should be higher than before", () => readingHigher(page, "tree zoom", el("grid")));
      await session.step(129, "When user scrolls the mouse wheel up 150 times over the \"tree\" area of grid holding Control", () => wheelOverAreaTimesHolding(page, "up", 150, "tree", el("grid"), "Control"));
      await session.step(130, "Then the \"tree zoom\" reading of grid should be 100", () => readingIs(page, "tree zoom", el("grid"), 100));
      await session.step(131, "When user scrolls the mouse wheel down 150 times over the \"tree\" area of grid holding Control", () => wheelOverAreaTimesHolding(page, "down", 150, "tree", el("grid"), "Control"));
      await session.step(132, "Then the \"tree zoom\" reading of grid should be 1", () => readingIs(page, "tree zoom", el("grid"), 1));
      await session.step(133, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reset Zoom from the tree's context menu puts the zoom back", async () => {
      await session.step(136, "When user scrolls the mouse wheel up 3 times over the \"tree\" area of grid holding Control", () => wheelOverAreaTimesHolding(page, "up", 3, "tree", el("grid"), "Control"));
      await session.step(137, "Then the \"tree zoom\" reading of grid should be higher than before", () => readingHigher(page, "tree zoom", el("grid")));
      await session.step(138, "When user picks \"Reset Zoom\" from the context menu of the \"tree\" area of grid", () => pickFromAreaContextMenu(page, "Reset Zoom", "tree", el("grid")));
      await session.step(139, "Then the \"tree zoom\" reading of grid should be 1", () => readingIs(page, "tree zoom", el("grid"), 1));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running the clustering again warns and replaces the tree", async () => {
      await session.step(143, "When user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
      await session.step(144, "And user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(145, "Then a warning balloon containing \"Closing existing dendrogram\" should have been shown", () => warningBalloonText(page, "Closing existing dendrogram"));
      await session.step(146, "And the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(147, "And there should be 1 visible \"Assign Clusters\" icon", () => visibleCount(page, 1, el("\"Assign Clusters\" icon")));
      await session.step(148, "And there should be 1 visible \"Remove Dendrogram\" icon", () => visibleCount(page, 1, el("\"Remove Dendrogram\" icon")));
      await session.step(149, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control with the wheel zooms the new tree in", async () => {
      await session.step(152, "When user scrolls the mouse wheel up 3 times over the \"tree\" area of grid holding Control", () => wheelOverAreaTimesHolding(page, "up", 3, "tree", el("grid"), "Control"));
      await session.step(153, "Then the \"tree zoom\" reading of grid should be higher than before", () => readingHigher(page, "tree zoom", el("grid")));
    });
    await run.scenario("A double-click on the tree's empty margin resets the zoom", async () => {
      await session.step(157, "When user double-clicks on the \"tree margin\" area of grid", () => doubleClickArea(page, "tree margin", el("grid")));
      await session.step(158, "Then the \"tree zoom\" reading of grid should be 1", () => readingIs(page, "tree zoom", el("grid"), 1));
    }, {knownFailure: true});
    run.finish();
  });
});
