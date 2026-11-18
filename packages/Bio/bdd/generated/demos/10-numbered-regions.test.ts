/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/10-numbered-regions.feature
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
import {carriesAtLeast} from '../../bindings/annotations.js';
import {gridColumnStaysWide} from '../../bindings/demo-antibody.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {balloonsClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noBalloons, noErrors, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridRowHeight, gridShowsColumns, headerStartsAt, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("6 · Antibodies and proteins · CDR3", () => {
  const session = feature(test, "features/demos/10-numbered-regions.feature", import.meta.url);
  test("CDR3 extracted with its IMGT positions", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(18, "And user opens antibodies dataset", () => openDataset(page, ds("antibodies")));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(23, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(24, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(25, "And the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(26, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(27, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned)\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned)"));
    await session.step(28, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"95\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "95"));
    await session.step(29, "And the \"AntibodyHC (aligned)\" column of the grid stays 900 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)", 900));
    await session.step(30, "Then the \"header positions of AntibodyHC (aligned)\" reading of grid should be \"100, 110, 120, C+10\"", () => readingReads(page, "header positions of AntibodyHC (aligned)", el("grid"), "100, 110, 120, C+10"));
    await session.step(31, "When user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 12 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 12 of AntibodyHC (aligned)", el("grid")));
    await session.step(32, "Then a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
    await session.step(33, "And the \"header tracks of AntibodyHC (aligned)(CDR3)\" reading of grid should be \"Conservation, WebLogo\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)(CDR3)", el("grid"), "Conservation, WebLogo"));
    await session.step(34, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(35, "And the rows of the grid are 30 px tall", () => gridRowHeight(page, 30));
    await session.step(36, "And the monomers of the \"AntibodyHC (aligned)(CDR3)\" column are drawn at 18 px", () => sequenceFontSize(page, "AntibodyHC (aligned)(CDR3)", 18));
    await session.step(37, "And the \"AntibodyHC (aligned)\" column of the grid stays 560 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)", 560));
    await session.step(38, "And the \"AntibodyHC (aligned)(CDR3)\" column of the grid stays 680 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)(CDR3)", 680));
    await session.step(39, "Then the \"header positions of AntibodyHC (aligned)(CDR3)\" reading of grid should be \"105, 110\"", () => readingReads(page, "header positions of AntibodyHC (aligned)(CDR3)", el("grid"), "105, 110"));
    await session.step(40, "When the video says \"CDR3 keeps its IMGT positions: 105 to 117\"", () => videoSays(page, "CDR3 keeps its IMGT positions: 105 to 117"));
    await session.step(41, "And user hovers over the \"position 111A of AntibodyHC (aligned)(CDR3) header\" area of grid", () => hoverArea(page, "position 111A of AntibodyHC (aligned)(CDR3) header", el("grid")));
    await session.step(42, "Then tooltip should contain text \"Position: 111A\"", () => shouldContainText(page, el("tooltip"), "Position: 111A"));
    await session.step(43, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."));
    await session.step(44, "And user selects \"AntibodyHC (aligned)(CDR3)\" in Sequence input in \"Scan Sequence Liabilities\" dialog", () => selectIn(page, "AntibodyHC (aligned)(CDR3)", el("Sequence input in \"Scan Sequence Liabilities\" dialog")));
    await session.step(45, "And user clicks on OK button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("OK button in \"Scan Sequence Liabilities\" dialog")));
    await session.step(46, "Then a new column \"~AntibodyHC (aligned)(CDR3)_annotations\" should have been added", () => newColumnNamed(page, "~AntibodyHC (aligned)(CDR3)_annotations"));
    await session.step(47, "And \"AntibodyHC (aligned)(CDR3)\" column should carry at least 3 annotations", () => carriesAtLeast(page, "AntibodyHC (aligned)(CDR3)", 3));
    await session.step(48, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(49, "When the video says \"Liabilities marked in every CDR3\"", () => videoSays(page, "Liabilities marked in every CDR3"));
    await session.step(50, "And user hovers over monomer 3 of the \"cell 12 of AntibodyHC (aligned)(CDR3)\" area of grid", () => hoverMonomer(page, 3, "cell 12 of AntibodyHC (aligned)(CDR3)", el("grid")));
    await session.step(51, "Then tooltip should contain text \"Isomerization (DG)\"", () => shouldContainText(page, el("tooltip"), "Isomerization (DG)"));
    await session.step(52, "When user hovers over the \"position 105 of AntibodyHC (aligned)(CDR3) header\" area of grid", () => hoverArea(page, "position 105 of AntibodyHC (aligned)(CDR3) header", el("grid")));
    await session.step(53, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
});
