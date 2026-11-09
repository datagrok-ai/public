/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/14-sar-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clusterStatistics} from '../../bindings/clusters.js';
import {headerStartsAt} from '../../bindings/demo-sar.js';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, selectIn, shouldBe, typeAtCaret, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {onlyOfSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {balloonsClosed, browsePanelClosed, browsePanelOpen, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, painted, pickFromAreaContextMenu, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns, viewersSettled} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Peptide SAR", () => {
  const session = feature(test, "features/demos/14-sar-viewers.feature", import.meta.url);
  test("Clusters, residues and cliffs in more viewers", {tag: ["@guide", "@demo", "@viewport:1920x1080"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(22, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---Bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio tree node inside browse tree")));
    await session.step(24, "And Files---App-Data---Bio---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio---samples tree node inside browse tree")));
    await session.step(25, "And user double-clicks on Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree")));
    await session.step(26, "And the \"antibodies\" view should be current", () => viewIsCurrent(page, "antibodies"));
    await session.step(27, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(28, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(29, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(30, "And user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(31, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(32, "And a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(33, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned), Y\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned), Y"));
    await session.step(34, "And the \"AntibodyHC (aligned)\" column of the grid is 1300 px wide", () => gridColumnWidth(page, "AntibodyHC (aligned)", 1300));
    await session.step(35, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"93\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "93"));
    await session.step(36, "And the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(37, "And user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 1 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 1 of AntibodyHC (aligned)", el("grid")));
    await session.step(38, "And a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
    await session.step(39, "And user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(40, "And user selects \"AntibodyHC (aligned)(CDR3)\" in Sequence input in \"Analyze Peptides\" dialog", () => selectIn(page, "AntibodyHC (aligned)(CDR3)", el("Sequence input in \"Analyze Peptides\" dialog")));
    await session.step(41, "And user selects \"-lg\" in Scaling input in \"Analyze Peptides\" dialog", () => selectIn(page, "-lg", el("Scaling input in \"Analyze Peptides\" dialog")));
    await session.step(42, "And user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(43, "And user enters \"90\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "90", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(44, "And user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(45, "And user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(46, "And the SAR analysis should be ready", () => sarReady(page));
    await session.step(47, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(48, "When the video says \"Each MCL cluster: its CDR3 logo and pKd\"", () => videoSays(page, "Each MCL cluster: its CDR3 logo and pKd"));
    await session.step(49, "And user hovers over the \"cluster 3\" area of Logo Summary Table viewer", () => hoverArea(page, "cluster 3", el("Logo Summary Table viewer")));
    await session.step(50, "Then the cluster statistics in tooltip should match cluster \"3\" of table \"antibodies\"", () => clusterStatistics(page, el("tooltip"), "3", "antibodies"));
    await session.step(51, "When user clicks on the \"cluster 3\" area of Logo Summary Table viewer", () => clickArea(page, "cluster 3", el("Logo Summary Table viewer")));
    await session.step(52, "Then only rows where \"Cluster (MCL)\" is \"3\" should be selected", () => onlyOfSelected(page, "Cluster (MCL)", "3"));
    await session.step(53, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(54, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(55, "When user types \"Position\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Position", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(56, "And user clicks on first \"Sequence Position Statistics\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sequence Position Statistics\" button in \"Add Viewer\" dialog")));
    await session.step(57, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(58, "And Sequence Position Statistics viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Position Statistics"));
    await session.step(59, "Given the viewers have settled", () => viewersSettled(page));
    await session.step(60, "When user clicks on \"Edit position\" icon in Sequence Position Statistics viewer", () => clickOn(page, el("\"Edit position\" icon in Sequence Position Statistics viewer")));
    await session.step(61, "And user types \"2\" at the caret", () => typeAtCaret(page, "2"));
    await session.step(62, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(63, "Then the \"positions\" reading of Sequence Position Statistics viewer should be \"2\"", () => readingReads(page, "positions", el("Sequence Position Statistics viewer"), "2"));
    await session.step(64, "When the video says \"pKd of every residue at one CDR3 position\"", () => videoSays(page, "pKd of every residue at one CDR3 position"));
    await session.step(65, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(66, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(67, "When user types \"Mutation\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Mutation", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(68, "And user clicks on first \"Sequence Mutation Cliffs\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sequence Mutation Cliffs\" button in \"Add Viewer\" dialog")));
    await session.step(69, "Then Sequence Mutation Cliffs viewer should be added to the open tableview", () => viewerAdded(page, "Sequence Mutation Cliffs"));
    await session.step(70, "When user selects \"8\" in Position input in Sequence Mutation Cliffs viewer", () => selectIn(page, "8", el("Position input in Sequence Mutation Cliffs viewer")));
    await session.step(71, "Then the \"position\" reading of Sequence Mutation Cliffs viewer should be 8", () => readingIs(page, "position", el("Sequence Mutation Cliffs viewer"), 8));
    await session.step(72, "And the \"cliff rows\" reading of Sequence Mutation Cliffs viewer should be 9", () => readingIs(page, "cliff rows", el("Sequence Mutation Cliffs viewer"), 9));
    await session.step(73, "And line chart viewer in Sequence Mutation Cliffs viewer should be painted", () => painted(page, el("line chart viewer in Sequence Mutation Cliffs viewer")));
    await session.step(74, "When the video says \"Cliffs at one position: W there costs affinity\"", () => videoSays(page, "Cliffs at one position: W there costs affinity"));
    await session.step(75, "Then no errors should have been logged", () => noErrors(page));
  });
});
