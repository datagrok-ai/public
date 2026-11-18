/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/18-clustering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {dockedRightOfGrid} from '../../bindings/demo-layout.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveValue, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newestMatchingDistinct, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical, colorLinearOver} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {contextPanelClosed, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickFirstArea, noBalloons, noErrors, painted, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("4 · A peptide series · Clustering", () => {
  const session = feature(test, "features/demos/18-clustering.feature", import.meta.url);
  test("Cluster the peptides into a tree and cut it", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(25, "And the \"sequence\" column of the grid is 270 px wide", () => gridColumnWidth(page, "sequence", 270));
    await session.step(26, "And the rows of the grid are 14 px tall", () => gridRowHeight(page, 14));
    await session.step(27, "And user colors \"activity\" column linearly from \"#4575b4\" to \"#d73027\" over 2.0 to 4.5", () => colorLinearOver(page, "activity", "#4575b4", "#d73027", 2, 4.5));
    await session.step(28, "And user watches the task bar", () => watchTaskBar(page));
    await session.step(29, "When user picks \"Bio > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Hierarchical Clustering..."));
    await session.step(30, "Then \"Hierarchical Clustering\" dialog should be visible", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "visible"));
    await session.step(31, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
    await session.step(32, "Then the task bar should have finished \"Creating dendrogram\"", () => taskBarFinished(page, "Creating dendrogram"));
    await session.step(33, "And the \"tree leaves\" reading of grid should be 99", () => readingIs(page, "tree leaves", el("grid"), 99));
    await session.step(34, "When the video says \"A tree of the peptides by sequence distance\"", () => videoSays(page, "A tree of the peptides by sequence distance"));
    await session.step(35, "And user clicks on the first \"node\" area of grid", () => clickFirstArea(page, "node", el("grid")));
    await session.step(36, "Then the \"tree selected leaves\" reading of grid should be at least 20", () => readingAtLeast(page, "tree selected leaves", el("grid"), 20));
    await session.step(37, "When user clicks on \"Assign Clusters\" icon", () => clickOn(page, el("\"Assign Clusters\" icon")));
    await session.step(38, "Then \"Assign Clusters\" dialog should be visible", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "visible"));
    await session.step(39, "When the video says \"Cut it where it splits into three\"", () => videoSays(page, "Cut it where it splits into three"));
    await session.step(40, "And user enters \"6.5\" into Threshold input in \"Assign Clusters\" dialog", () => enterInto(page, "6.5", el("Threshold input in \"Assign Clusters\" dialog")));
    await session.step(41, "Then Clusters input in \"Assign Clusters\" dialog should have value \"3\"", () => shouldHaveValue(page, el("Clusters input in \"Assign Clusters\" dialog"), "3"));
    await session.step(42, "When user unchecks \"Medoid columns\" checkbox in \"Assign Clusters\" dialog", () => uncheck(page, el("\"Medoid columns\" checkbox in \"Assign Clusters\" dialog")));
    await session.step(43, "And user clicks on Assign button in \"Assign Clusters\" dialog", () => clickOn(page, el("Assign button in \"Assign Clusters\" dialog")));
    await session.step(44, "Then the \"Assign Clusters\" dialog should close", () => dialogCloses(page, "Assign Clusters"));
    await session.step(45, "And the newest column matching \"^Cluster \\(\" should have 3 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(", 3));
    await session.step(46, "And every value of \"Cluster (6.50)\" column should match \"^[123]$\"", () => everyValueMatches(page, "Cluster (6.50)", "^[123]$"));
    await session.step(47, "And the table should not have a column \"Medoid Rank (6.50)\"", () => hasNoColumn(page, "Medoid Rank (6.50)"));
    await session.step(48, "Given the grid shows only the columns \"sequence, activity, Cluster (6.50)\"", () => gridShowsColumns(page, "sequence, activity, Cluster (6.50)"));
    await session.step(49, "And user colors \"Cluster (6.50)\" column categorically:", () => colorCategorical(page, "Cluster (6.50)", [["1","#8e6cc0"],["2","#2a9d8f"],["3","#e9c46a"]]), [["1","#8e6cc0"],["2","#2a9d8f"],["3","#e9c46a"]]);
    await session.step(53, "And the \"Cluster (6.50)\" column of the grid is 125 px wide", () => gridColumnWidth(page, "Cluster (6.50)", 125));
    await session.step(54, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["valueColumnName","activity"],["category1ColumnName","Cluster (6.50)"]]), [["valueColumnName","activity"],["category1ColumnName","Cluster (6.50)"]]);
    await session.step(57, "And box plot viewer is docked to the right of the grid, 45% of the view wide", () => dockedRightOfGrid(page, el("box plot viewer"), 45));
    await session.step(58, "And user scrolls the grid to the \"sequence\" column", () => scrollGridTo(page, "sequence"));
    await session.step(59, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(60, "When the video says \"Three clusters, three families: one stands out in activity\"", () => videoSays(page, "Three clusters, three families: one stands out in activity"));
    await session.step(61, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
});
