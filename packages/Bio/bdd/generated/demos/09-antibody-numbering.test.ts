/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/09-antibody-numbering.feature
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
import {dragRuler, gridColumnStaysWide, hoverRegionTrack} from '../../bindings/demo-antibody.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {balloonsClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden, tooltipHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns, headerStartsAt, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("6 · Antibodies and proteins · Numbering", () => {
  const session = feature(test, "features/demos/09-antibody-numbering.feature", import.meta.url);
  test("Antibody numbering in the browser", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(20, "And user opens antibodies dataset", () => openDataset(page, ds("antibodies")));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the grid shows only the columns \"Antibody_ID, AntibodyHC\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC"));
    await session.step(25, "And the \"AntibodyHC\" column of the grid is 1400 px wide", () => gridColumnWidth(page, "AntibodyHC", 1400));
    await session.step(26, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(27, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
    await session.step(28, "And Engine input in \"Apply Antibody Numbering\" dialog should have value \"ANARCI (HMMER)\"", () => shouldHaveValue(page, el("Engine input in \"Apply Antibody Numbering\" dialog"), "ANARCI (HMMER)"));
    await session.step(29, "And Scheme input in \"Apply Antibody Numbering\" dialog should have value \"imgt\"", () => shouldHaveValue(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt"));
    await session.step(30, "And Scheme input in \"Apply Antibody Numbering\" dialog should offer \"imgt, kabat, chothia, martin, aho, wolfguy\"", () => shouldOffer(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt, kabat, chothia, martin, aho, wolfguy"));
    await session.step(31, "When user selects \"imgt\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "imgt", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
    await session.step(32, "Then Scheme input in \"Apply Antibody Numbering\" dialog should have value \"imgt\"", () => shouldHaveValue(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt"));
    await session.step(33, "When user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(34, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(35, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "imgt"));
    await session.step(36, "And \"AntibodyHC\" column should carry at least 7 annotations", () => carriesAtLeast(page, "AntibodyHC", 7));
    await session.step(37, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(38, "And the \"AntibodyHC\" column of the grid is 1400 px wide", () => gridColumnWidth(page, "AntibodyHC", 1400));
    await session.step(39, "When the video says \"Frameworks and CDRs, colored on every chain\"", () => videoSays(page, "Frameworks and CDRs, colored on every chain"));
    await session.step(40, "And user hovers over monomer 54 of the \"cell 3 of AntibodyHC\" area of grid", () => hoverMonomer(page, 54, "cell 3 of AntibodyHC", el("grid")));
    await session.step(41, "Then tooltip should contain text \"CDR2\"", () => shouldContainText(page, el("tooltip"), "CDR2"));
    await session.step(42, "When the video says \"The aligned column, ruled in IMGT positions\"", () => videoSays(page, "The aligned column, ruled in IMGT positions"));
    await session.step(43, "Given the tooltip is hidden", () => tooltipHidden(page));
    await session.step(44, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned)\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned)"));
    await session.step(45, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"1\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "1"));
    await session.step(46, "Then the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(47, "Given the \"AntibodyHC (aligned)\" column of the grid stays 1400 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)", 1400));
    await session.step(48, "Then the \"header positions of AntibodyHC (aligned)\" reading of grid should be \"10, 20, 30, 40, 50, 60, 70, 80, 90\"", () => readingReads(page, "header positions of AntibodyHC (aligned)", el("grid"), "10, 20, 30, 40, 50, 60, 70, 80, 90"));
    await session.step(49, "When user drags the ruler under the \"position 85 of AntibodyHC (aligned) header\" area of grid by 1100 pixels to the left", () => dragRuler(page, "position 85 of AntibodyHC (aligned) header", el("grid"), 1100, "left"));
    await session.step(50, "Then the \"header positions of AntibodyHC (aligned)\" reading of grid should be \"80, 90, 100, 110, 120, C+10, C+20\"", () => readingReads(page, "header positions of AntibodyHC (aligned)", el("grid"), "80, 90, 100, 110, 120, C+10, C+20"));
    await session.step(51, "When user hovers over the \"position 110 of AntibodyHC (aligned) header\" area of grid", () => hoverArea(page, "position 110 of AntibodyHC (aligned) header", el("grid")));
    await session.step(52, "Then tooltip should contain text \"Position: 110\"", () => shouldContainText(page, el("tooltip"), "Position: 110"));
    await session.step(53, "When user hovers over the region track under the \"position 111A of AntibodyHC (aligned) header\" area of grid", () => hoverRegionTrack(page, "position 111A of AntibodyHC (aligned) header", el("grid")));
    await session.step(54, "Then tooltip should contain text \"CDR3\"", () => shouldContainText(page, el("tooltip"), "CDR3"));
    await session.step(55, "And tooltip should contain text \"Position: 111A\"", () => shouldContainText(page, el("tooltip"), "Position: 111A"));
    await session.step(56, "When the video says \"Light chains, numbered the same way\"", () => videoSays(page, "Light chains, numbered the same way"));
    await session.step(57, "Given user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(58, "And user selects \"AntibodyLC\" in Sequence input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "AntibodyLC", el("Sequence input in \"Apply Antibody Numbering\" dialog")));
    await session.step(59, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(60, "Then a new column \"AntibodyLC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyLC (aligned)"));
    await session.step(61, "And \"AntibodyLC (aligned)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyLC (aligned)", ".numberingScheme", "imgt"));
    await session.step(62, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(63, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned), AntibodyLC (aligned)\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned), AntibodyLC (aligned)"));
    await session.step(64, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"100\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "100"));
    await session.step(65, "And the header of the \"AntibodyLC (aligned)\" column starts at position \"100\"", () => headerStartsAt(page, "AntibodyLC (aligned)", "100"));
    await session.step(66, "Then the \"header tracks of AntibodyLC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyLC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(67, "Given the \"AntibodyHC (aligned)\" column of the grid stays 700 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)", 700));
    await session.step(68, "And the \"AntibodyLC (aligned)\" column of the grid stays 700 px wide", () => gridColumnStaysWide(page, "AntibodyLC (aligned)", 700));
    await session.step(69, "Then the \"header positions of AntibodyLC (aligned)\" reading of grid should be \"100, 110, 120\"", () => readingReads(page, "header positions of AntibodyLC (aligned)", el("grid"), "100, 110, 120"));
    await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
});
