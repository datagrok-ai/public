/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/11b-profile-hmm.feature
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
import {gridColumnStaysWide, headerStartsAt} from '../../bindings/demo-antibody.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectFirstRows, selectNoRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {balloonsClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Profile HMMs from your data", () => {
  const session = feature(test, "features/demos/11b-profile-hmm.feature", import.meta.url);
  test("A profile HMM built from selected chains", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(16, "And user opens antibodies dataset", () => openDataset(page, ds("antibodies")));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(21, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(22, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(23, "And the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
    await session.step(24, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(25, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned)\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned)"));
    await session.step(26, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"1\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "1"));
    await session.step(27, "And the \"AntibodyHC (aligned)\" column of the grid stays 1400 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)", 1400));
    await session.step(28, "When user selects the first 15 rows", () => selectFirstRows(page, 15));
    await session.step(29, "And the video says \"15 chains picked to build a profile HMM\"", () => videoSays(page, "15 chains picked to build a profile HMM"));
    await session.step(30, "And user picks \"Bio > Search > Build Profile HMM...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Build Profile HMM..."));
    await session.step(31, "Then \"Build Profile HMM\" dialog should be visible", () => shouldBe(page, el("\"Build Profile HMM\" dialog"), "visible"));
    await session.step(32, "And Rows input in \"Build Profile HMM\" dialog should have value \"Selected rows\"", () => shouldHaveValue(page, el("Rows input in \"Build Profile HMM\" dialog"), "Selected rows"));
    await session.step(33, "When user enters \"my15\" into Name input in \"Build Profile HMM\" dialog", () => enterInto(page, "my15", el("Name input in \"Build Profile HMM\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Build Profile HMM\" dialog", () => clickOn(page, el("OK button in \"Build Profile HMM\" dialog")));
    await session.step(35, "Then a new column \"my15 E-value\" should have been added", () => newColumnNamed(page, "my15 E-value"));
    await session.step(36, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(37, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(38, "And the grid shows only the columns \"Antibody_ID, AntibodyHC (aligned), my15 E-value\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC (aligned), my15 E-value"));
    await session.step(39, "And the header of the \"AntibodyHC (aligned)\" column starts at position \"1\"", () => headerStartsAt(page, "AntibodyHC (aligned)", "1"));
    await session.step(40, "And the \"AntibodyHC (aligned)\" column of the grid stays 1150 px wide", () => gridColumnStaysWide(page, "AntibodyHC (aligned)", 1150));
    await session.step(41, "And the \"my15 E-value\" column of the grid is 130 px wide", () => gridColumnWidth(page, "my15 E-value", 130));
    await session.step(42, "When the video says \"The profile scores every chain; its logo selects them\"", () => videoSays(page, "The profile scores every chain; its logo selects them"));
    await session.step(43, "And user selects no rows", () => selectNoRows(page));
    await session.step(44, "And user hovers over the \"monomer R at position 38\" area of WebLogo viewer", () => hoverArea(page, "monomer R at position 38", el("WebLogo viewer")));
    await session.step(45, "And user clicks on the \"monomer R at position 38\" area of WebLogo viewer", () => clickArea(page, "monomer R at position 38", el("WebLogo viewer")));
    await session.step(46, "Then 378 rows should be selected", () => selectedRowCount(page, 378));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
