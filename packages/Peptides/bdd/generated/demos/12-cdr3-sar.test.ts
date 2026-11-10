/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/12-cdr3-sar.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {headerStartsAt} from '../../bindings/demo-sar.js';
import {peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelClosed, browsePanelOpen, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, painted, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Peptide SAR", () => {
  const session = feature(test, "features/demos/12-cdr3-sar.feature", import.meta.url);
  test("From antibodies to a CDR3 SAR", {tag: ["@guide", "@demo", "@viewport:1920x1080"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(19, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(20, "And Files---App-Data---Bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio tree node inside browse tree")));
    await session.step(21, "And Files---App-Data---Bio---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio---samples tree node inside browse tree")));
    await session.step(22, "When user double-clicks on Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree")));
    await session.step(23, "Then the \"antibodies\" view should be current", () => viewIsCurrent(page, "antibodies"));
    await session.step(24, "And \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await session.step(25, "Given the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(27, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(28, "And the grid shows only the columns \"Antibody_ID, AntibodyHC, Y\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC, Y"));
    await session.step(29, "And the \"AntibodyHC\" column of the grid is 1300 px wide", () => gridColumnWidth(page, "AntibodyHC", 1300));
    await session.step(30, "And user sets the format of \"Y\" column in the grid to \"scientific\"", () => setGridColumnFormat(page, "Y", "scientific"));
    await session.step(31, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(32, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
    await session.step(33, "When user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(34, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(35, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "imgt"));
    await session.step(36, "Given the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned), Y\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned), Y"));
    await session.step(37, "And the \"AntibodyHC (aligned)\" column of the grid is 1300 px wide", () => gridColumnWidth(page, "AntibodyHC (aligned)", 1300));
    await session.step(38, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"93\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "93"));
    await session.step(39, "Then the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(40, "When the video says \"CDR3, the loop that binds the antigen, becomes a column\"", () => videoSays(page, "CDR3, the loop that binds the antigen, becomes a column"));
    await session.step(41, "And user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 1 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 1 of AntibodyHC (aligned)", el("grid")));
    await session.step(42, "Then a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
    await session.step(43, "Given the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned), AntibodyHC (aligned)(CDR3), Y\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned), AntibodyHC (aligned)(CDR3), Y"));
    await session.step(44, "And the \"AntibodyHC (aligned)\" column of the grid is 820 px wide", () => gridColumnWidth(page, "AntibodyHC (aligned)", 820));
    await session.step(45, "And the \"AntibodyHC (aligned)(CDR3)\" column of the grid is 560 px wide", () => gridColumnWidth(page, "AntibodyHC (aligned)(CDR3)", 560));
    await session.step(46, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(47, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(48, "When user selects \"AntibodyHC (aligned)(CDR3)\" in Sequence input in \"Analyze Peptides\" dialog", () => selectIn(page, "AntibodyHC (aligned)(CDR3)", el("Sequence input in \"Analyze Peptides\" dialog")));
    await session.step(49, "And the video says \"Y is a Kd in molar: -lg turns it into pKd\"", () => videoSays(page, "Y is a Kd in molar: -lg turns it into pKd"));
    await session.step(50, "And user selects \"-lg\" in Scaling input in \"Analyze Peptides\" dialog", () => selectIn(page, "-lg", el("Scaling input in \"Analyze Peptides\" dialog")));
    await session.step(51, "And user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(52, "And user enters \"90\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "90", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(53, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(54, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(55, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(56, "And the SAR setting \"sequenceColumnName\" should be \"AntibodyHC (aligned)(CDR3)\"", () => sarSetting(page, "sequenceColumnName", "AntibodyHC (aligned)(CDR3)"));
    await session.step(57, "And the SAR setting \"activityScaling\" should be \"-lg\"", () => sarSetting(page, "activityScaling", "-lg"));
    await session.step(58, "And the \"completed threshold\" reading of MCL viewer should be 90", () => readingIs(page, "completed threshold", el("MCL viewer"), 90));
    await session.step(59, "And the \"positions\" reading of Sequence Variability Map viewer should be 26", () => readingIs(page, "positions", el("Sequence Variability Map viewer"), 26));
    await session.step(60, "And the \"clusters shown\" reading of Logo Summary Table viewer should be 11", () => readingIs(page, "clusters shown", el("Logo Summary Table viewer"), 11));
    await session.step(61, "And Sequence Variability Map viewer should be painted", () => painted(page, el("Sequence Variability Map viewer")));
    await session.step(62, "And Most Potent Residues viewer should be painted", () => painted(page, el("Most Potent Residues viewer")));
    await session.step(63, "And Logo Summary Table viewer should be painted", () => painted(page, el("Logo Summary Table viewer")));
    await session.step(64, "When the video says \"Variability map, most potent residues, MCL clusters and their CDR3 logos\"", () => videoSays(page, "Variability map, most potent residues, MCL clusters and their CDR3 logos"));
    await session.step(65, "Then no errors should have been logged", () => noErrors(page));
  });
});
