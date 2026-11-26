/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/08-hierarchical-clustering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldContainText, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnMatching, newestMatchingDistinct, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical, colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelClosed, contextPanelClosed, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, hoverArea, noBalloons, noErrors, pickFromAreaContextMenu, pointerAway, readingFinite, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridRowHeight, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/08-hierarchical-clustering.feature", import.meta.url);
  test("Cluster rows into a tree", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And simple mode is off", () => simpleModeOff(page));
    await session.step(26, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(27, "And user colors \"Species\" column categorically:", () => colorCategorical(page, "Species", [["setosa","#1f77b4"],["versicolor","#2ca02c"],["virginica","#ffbb78"]]), [["setosa","#1f77b4"],["versicolor","#2ca02c"],["virginica","#ffbb78"]]);
    await session.step(31, "And the rows of the grid are 6 px tall", () => gridRowHeight(page, 6));
    await session.step(32, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(33, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(34, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(35, "When the video says \"Group the flowers by their four measurements\" at the bottom", () => videoSaysAt(page, "Group the flowers by their four measurements", "bottom"));
    await session.step(37, "And user picks \"ML > Cluster > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > Hierarchical Clustering..."), undefined, "Hierarchical Clustering is under ML > Cluster");
    await session.step(38, "Then \"Hierarchical Clustering\" dialog should be visible", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "visible"));
    await session.step(40, "When user clicks on Features input in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("Features input in \"Hierarchical Clustering\" dialog")), undefined, "Features: the columns that make rows alike");
    await session.step(41, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(42, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
    await session.step(44, "And user toggles the \"col 1\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "col 1", el("\"Select columns...\" dialog")), undefined, "Leave out the row number");
    await session.step(45, "Then \"Select columns...\" dialog should contain text \"4 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "4 checked"));
    await session.step(46, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(48, "And user hovers over Linkage input in \"Hierarchical Clustering\" dialog", () => hoverOver(page, el("Linkage input in \"Hierarchical Clustering\" dialog")), undefined, "Distance and Linkage: how closeness is measured and groups are joined");
    await session.step(49, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(50, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
    await session.step(51, "Then the task bar should have finished \"Creating dendrogram\"", () => taskBarFinished(page, "Creating dendrogram"));
    await session.step(52, "And the \"tree leaves\" reading of grid should be 150", () => readingIs(page, "tree leaves", el("grid"), 150));
    await session.step(53, "When the video says \"The tree beside the grid: one leaf per row, branches follow the species\" at the bottom", () => videoSaysAt(page, "The tree beside the grid: one leaf per row, branches follow the species", "bottom"));
    await session.step(54, "And user hovers over the \"tree\" area of grid", () => hoverArea(page, "tree", el("grid")));
    await session.step(56, "And user picks \"Assign Clusters\" from the context menu of the \"tree\" area of grid", () => pickFromAreaContextMenu(page, "Assign Clusters", "tree", el("grid")), undefined, "Right-click the tree to cut it into clusters");
    await session.step(57, "Then \"Assign Clusters\" dialog should be visible", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "visible"));
    await session.step(59, "When user enters \"3\" into Clusters input in \"Assign Clusters\" dialog", () => enterInto(page, "3", el("Clusters input in \"Assign Clusters\" dialog")), undefined, "Clusters: how many groups to cut the tree into");
    await session.step(60, "Then the \"cut threshold\" reading of grid should be a finite number", () => readingFinite(page, "cut threshold", el("grid")));
    await session.step(61, "When user clicks on Assign button in \"Assign Clusters\" dialog", () => clickOn(page, el("Assign button in \"Assign Clusters\" dialog")));
    await session.step(62, "Then the \"Assign Clusters\" dialog should close", () => dialogCloses(page, "Assign Clusters"));
    await session.step(63, "And a new column matching \"^Cluster \\(\" should have been added", () => newColumnMatching(page, "^Cluster \\("));
    await session.step(64, "And the newest column matching \"^Cluster \\(\" should have 3 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(", 3));
    await session.step(66, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header Cluster (38.58)\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header Cluster (38.58)", el("grid")), undefined, "Assign adds a Cluster column; color it to see the groups");
    await session.step(67, "Then \"Cluster (38.58)\" column should be color-coded categorically", () => colorCodedCategorically(page, "Cluster (38.58)"));
    await session.step(68, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(70, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")), undefined, "Setosa gets a cluster of its own; the other two species overlap");
    await session.step(71, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(72, "And no errors should have been logged", () => noErrors(page));
  });
});
