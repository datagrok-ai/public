/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/16b-chemical-space.feature
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
import {dockedRightOfViewer} from '../../bindings/demo-layout.js';
import {boxAroundMarkers} from '../../bindings/demo-space.js';
import {atomicLevelOffCamera, bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, newColumnMatching, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical, noneOfSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("4 · A peptide series · Chemical space", () => {
  const session = feature(test, "features/demos/16b-chemical-space.feature", import.meta.url);
  test("The same peptides in chemical space", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(23, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(24, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(25, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(26, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(27, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(28, "And the grid shows only the columns \"sequence\"", () => gridShowsColumns(page, "sequence"));
    await session.step(29, "And the \"sequence\" column of the grid is 250 px wide", () => gridColumnWidth(page, "sequence", 250));
    await session.step(30, "And user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
    await session.step(31, "And user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(32, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(33, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(34, "Given the \"sequence\" column is converted to atomic level drawn straight with its monomers highlighted", () => atomicLevelOffCamera(page, "sequence"));
    await session.step(35, "And the grid shows only the columns \"sequence\"", () => gridShowsColumns(page, "sequence"));
    await session.step(36, "When user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."));
    await session.step(37, "Then \"Chem Space\" dialog should be visible", () => shouldBe(page, el("\"Chem Space\" dialog"), "visible"));
    await session.step(38, "And Column input in \"Chem Space\" dialog should contain text \"molfile(sequence)\"", () => shouldContainText(page, el("Column input in \"Chem Space\" dialog"), "molfile(sequence)"));
    await session.step(39, "When user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")));
    await session.step(40, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(41, "And a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
    await session.step(42, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
    await session.step(43, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
    await session.step(44, "Given second scatter plot viewer is docked to the right of scatter plot viewer", () => dockedRightOfViewer(page, el("second scatter plot viewer"), el("scatter plot viewer")));
    await session.step(45, "And the grid shows only the columns \"sequence\"", () => gridShowsColumns(page, "sequence"));
    await session.step(46, "And the rows of the grid are 25 px tall", () => gridRowHeight(page, 25));
    await session.step(47, "And user colors \"Cluster (DBSCAN)\" column categorically:", () => colorCategorical(page, "Cluster (DBSCAN)", [["1","#8e6cc0"],["2","#2a9d8f"],["3","#4575b4"]]), [["1","#8e6cc0"],["2","#2a9d8f"],["3","#4575b4"]]);
    await session.step(51, "And user sets \"colorColumnName\" property of second scatter plot viewer to \"Cluster (DBSCAN)\"", () => setProperty(page, "colorColumnName", el("second scatter plot viewer"), "Cluster (DBSCAN)"));
    await session.step(52, "And user sets \"displayLabels\" property of second scatter plot viewer to \"Never\"", () => setProperty(page, "displayLabels", el("second scatter plot viewer"), "Never"));
    await session.step(53, "Then second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
    await session.step(54, "When the video says \"Selected by sequence, the same peptides light up in chemical space\"", () => videoSays(page, "Selected by sequence, the same peptides light up in chemical space"));
    await session.step(55, "And user drags a selection box around the markers of rows 1 to 20 of scatter plot viewer", () => boxAroundMarkers(page, 1, 20, el("scatter plot viewer")));
    await session.step(56, "Then some rows should be selected", () => someSelected(page));
    await session.step(57, "And no rows where \"cluster\" is \"1\" should be selected", () => noneOfSelected(page, "cluster", "1"));
    await session.step(58, "And no rows where \"cluster\" is \"2\" should be selected", () => noneOfSelected(page, "cluster", "2"));
    await session.step(59, "And second scatter plot viewer should show a selection highlight", () => someHighlight(page, el("second scatter plot viewer")));
    await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
  });
});
