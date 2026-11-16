/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/13-mutation-cliffs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {balloonsClosed, browsePanelClosed, browsePanelOpen, contextPanelOpen, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge, dragTopBorder, gridColumnWidth, gridShowsColumns, headerStartsAt, viewersSettled} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("6 · Antibodies and proteins · Mutation cliffs", () => {
  const session = feature(test, "features/demos/13-mutation-cliffs.feature", import.meta.url);
  test("Mutation cliffs and the invariant map", {tag: ["@guide", "@demo", "@viewport:1920x1080"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And simple mode is off", () => simpleModeOff(page));
    await session.step(23, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(27, "And Files---App-Data---Bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio tree node inside browse tree")));
    await session.step(28, "And Files---App-Data---Bio---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio---samples tree node inside browse tree")));
    await session.step(29, "And user double-clicks on Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree")));
    await session.step(30, "And the \"antibodies\" view should be current", () => viewIsCurrent(page, "antibodies"));
    await session.step(31, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(32, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(33, "And user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(34, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(35, "And a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(36, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned), Y\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned), Y"));
    await session.step(37, "And the \"AntibodyHC (aligned)\" column of the grid is 1300 px wide", () => gridColumnWidth(page, "AntibodyHC (aligned)", 1300));
    await session.step(38, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"93\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "93"));
    await session.step(39, "And the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(40, "And user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 1 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 1 of AntibodyHC (aligned)", el("grid")));
    await session.step(41, "And a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
    await session.step(42, "And user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(43, "And user selects \"AntibodyHC (aligned)(CDR3)\" in Sequence input in \"Analyze Peptides\" dialog", () => selectIn(page, "AntibodyHC (aligned)(CDR3)", el("Sequence input in \"Analyze Peptides\" dialog")));
    await session.step(44, "And user selects \"-lg\" in Scaling input in \"Analyze Peptides\" dialog", () => selectIn(page, "-lg", el("Scaling input in \"Analyze Peptides\" dialog")));
    await session.step(45, "And user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(46, "And user enters \"90\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "90", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(47, "And user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(48, "And user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(49, "And the SAR analysis should be ready", () => sarReady(page));
    await session.step(50, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(51, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(52, "And user docks Sequence Variability Map viewer to the bottom edge of the view", () => dockToViewEdge(page, el("Sequence Variability Map viewer"), "bottom"));
    await session.step(53, "And user drags the top border of Sequence Variability Map viewer by 200 pixels up", () => dragTopBorder(page, el("Sequence Variability Map viewer"), 200));
    await session.step(54, "And the viewers have settled", () => viewersSettled(page));
    await session.step(55, "When user hovers over the \"cell W at 8\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell W at 8", el("Sequence Variability Map viewer")));
    await session.step(56, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(57, "And \"Pairs count\" table row in tooltip should have text \"Pairs count13\"", () => shouldHaveText(page, el("\"Pairs count\" table row in tooltip"), "Pairs count13"));
    await session.step(58, "And \"Mean activity\" table row in tooltip should not contain text \"0.000\"", () => shouldNotContainText(page, el("\"Mean activity\" table row in tooltip"), "0.000"));
    await session.step(59, "When the video says \"Mutation cliffs: CDR3 pairs one substitution apart\"", () => videoSays(page, "Mutation cliffs: CDR3 pairs one substitution apart"));
    await session.step(60, "And user clicks on the \"cell W at 8\" area of Sequence Variability Map viewer", () => clickArea(page, "cell W at 8", el("Sequence Variability Map viewer")));
    await session.step(61, "Then the \"selected monomer-positions\" reading of Sequence Variability Map viewer should be \"8:W\"", () => readingReads(page, "selected monomer-positions", el("Sequence Variability Map viewer"), "8:W"));
    await session.step(62, "And context panel should contain text \"Mutation Cliffs pairs\"", () => shouldContainText(page, el("context panel"), "Mutation Cliffs pairs"));
    await session.step(63, "When user hovers over \"Mutation Cliffs pairs\" pane in context panel", () => hoverOver(page, el("\"Mutation Cliffs pairs\" pane in context panel")));
    await session.step(64, "And user clicks on \"Expand to full screen\" icon in context panel", () => clickOn(page, el("\"Expand to full screen\" icon in context panel")));
    await session.step(65, "Then \"Mutation Cliffs pairs\" dialog should be visible", () => shouldBe(page, el("\"Mutation Cliffs pairs\" dialog"), "visible"));
    await session.step(66, "When the video says \"Every pair: the substitution and its jump in pKd\"", () => videoSays(page, "Every pair: the substitution and its jump in pKd"));
    await session.step(67, "And user hovers over the \"cell 6 of Mutation\" area of grid in \"Mutation Cliffs pairs\" dialog", () => hoverArea(page, "cell 6 of Mutation", el("grid in \"Mutation Cliffs pairs\" dialog")));
    await session.step(68, "And user hovers over the \"cell 6 of Delta\" area of grid in \"Mutation Cliffs pairs\" dialog", () => hoverArea(page, "cell 6 of Delta", el("grid in \"Mutation Cliffs pairs\" dialog")));
    await session.step(69, "And user clicks on CANCEL button in \"Mutation Cliffs pairs\" dialog", () => clickOn(page, el("CANCEL button in \"Mutation Cliffs pairs\" dialog")));
    await session.step(70, "Then \"Mutation Cliffs pairs\" dialog should be absent", () => shouldBe(page, el("\"Mutation Cliffs pairs\" dialog"), "absent"));
    await session.step(71, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
    await session.step(72, "Then the \"mode\" reading of Sequence Variability Map viewer should be \"Invariant Map\"", () => readingReads(page, "mode", el("Sequence Variability Map viewer"), "Invariant Map"));
    await session.step(73, "When user hovers over the \"cell W at 6\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell W at 6", el("Sequence Variability Map viewer")));
    await session.step(74, "Then Count table row in tooltip should contain text \"6 (\"", () => shouldContainText(page, el("Count table row in tooltip"), "6 ("));
    await session.step(75, "When the video says \"Invariant map: residues colored by mean pKd; a click selects their antibodies\"", () => videoSays(page, "Invariant map: residues colored by mean pKd; a click selects their antibodies"));
    await session.step(76, "And user clicks on the \"cell W at 6\" area of Sequence Variability Map viewer", () => clickArea(page, "cell W at 6", el("Sequence Variability Map viewer")));
    await session.step(77, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
    await session.step(78, "And context panel should contain text \"Invariant map: 6:W\"", () => shouldContainText(page, el("context panel"), "Invariant map: 6:W"));
    await session.step(79, "And context panel should contain text \"Mutation cliffs: 8:W\"", () => shouldContainText(page, el("context panel"), "Mutation cliffs: 8:W"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
});
