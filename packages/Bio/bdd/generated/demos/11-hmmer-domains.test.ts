/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/11-hmmer-domains.feature
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
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, selectIn, shouldBe, shouldContainText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {someValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {balloonsClosed, contextPanelClosed, contextPanelOpen, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, propertyShouldBe, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns, headerStartsAt} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("6 · Antibodies and proteins · Pfam domains", () => {
  const session = feature(test, "features/demos/11-hmmer-domains.feature", import.meta.url);
  test("Pfam domains, germlines and profile logos", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens antibodies dataset", () => openDataset(page, ds("antibodies")));
    await session.step(20, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the grid shows only the columns \"Antibody_ID, AntibodyHC\"", () => gridShowsColumns(page, "Antibody_ID, AntibodyHC"));
    await session.step(24, "And the \"AntibodyHC\" column of the grid is 1400 px wide", () => gridColumnWidth(page, "AntibodyHC", 1400));
    await session.step(25, "And the header of the \"AntibodyHC\" column starts at position \"85\"", () => headerStartsAt(page, "AntibodyHC", "85"));
    await session.step(26, "When user picks \"Bio > Annotate > Find Domains (HMMER)...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Find Domains (HMMER)..."));
    await session.step(27, "Then \"Find Domains (HMMER)\" dialog should be visible", () => shouldBe(page, el("\"Find Domains (HMMER)\" dialog"), "visible"));
    await session.step(28, "When user clicks on OK button in \"Find Domains (HMMER)\" dialog", () => clickOn(page, el("OK button in \"Find Domains (HMMER)\" dialog")));
    await session.step(29, "Then a new column \"AntibodyHC domains\" should have been added", () => newColumnNamed(page, "AntibodyHC domains"));
    await session.step(30, "And some value of \"AntibodyHC domains\" column should contain \"V-set\"", () => someValueContains(page, "AntibodyHC domains", "V-set"));
    await session.step(31, "And some value of \"AntibodyHC domains\" column should contain \"C1-set\"", () => someValueContains(page, "AntibodyHC domains", "C1-set"));
    await session.step(32, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(33, "And the \"AntibodyHC\" column of the grid is 1100 px wide", () => gridColumnWidth(page, "AntibodyHC", 1100));
    await session.step(34, "And the \"AntibodyHC domains\" column of the grid is 300 px wide", () => gridColumnWidth(page, "AntibodyHC domains", 300));
    await session.step(35, "When the video says \"Pfam domains of every chain: V-set, then C1-set\"", () => videoSays(page, "Pfam domains of every chain: V-set, then C1-set"));
    await session.step(36, "And user clicks on the \"cell 3 of AntibodyHC\" area of grid", () => clickArea(page, "cell 3 of AntibodyHC", el("grid")));
    await session.step(37, "Then the \"current row\" reading of grid should be 3", () => readingIs(page, "current row", el("grid"), 3));
    await session.step(38, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(39, "And the \"AntibodyHC\" column of the grid is 700 px wide", () => gridColumnWidth(page, "AntibodyHC", 700));
    await session.step(40, "And the header of the \"AntibodyHC\" column starts at position \"100\"", () => headerStartsAt(page, "AntibodyHC", "100"));
    await session.step(41, "When user expands Bioinformatics accordion header in context panel", () => expand(page, el("Bioinformatics accordion header in context panel")));
    await session.step(42, "And user expands ANARCI accordion header in context panel", () => expand(page, el("ANARCI accordion header in context panel")));
    await session.step(43, "Then \"ANARCI\" pane in context panel should contain the text \"IGHV\"", () => shouldContainText(page, el("\"ANARCI\" pane in context panel"), "IGHV"));
    await session.step(44, "When the video says \"ANARCI and Pfam panes: germlines, CDRs and domains\"", () => videoSays(page, "ANARCI and Pfam panes: germlines, CDRs and domains"));
    await session.step(45, "When user expands \"Pfam Domains\" accordion header in context panel", () => expand(page, el("\"Pfam Domains\" accordion header in context panel")));
    await session.step(46, "Then \"Pfam Domains\" pane in context panel should contain the text \"V-set\"", () => shouldContainText(page, el("\"Pfam Domains\" pane in context panel"), "V-set"));
    await session.step(47, "When user picks \"Bio > Analyze > Profile HMM Logo...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Profile HMM Logo..."));
    await session.step(48, "Then \"Profile HMM Logo\" dialog should be visible", () => shouldBe(page, el("\"Profile HMM Logo\" dialog"), "visible"));
    await session.step(49, "When user selects \"V-set\" in Family input in \"Profile HMM Logo\" dialog", () => selectIn(page, "V-set", el("Family input in \"Profile HMM Logo\" dialog")));
    await session.step(50, "Then Family input in \"Profile HMM Logo\" dialog should have value \"V-set\"", () => shouldHaveValue(page, el("Family input in \"Profile HMM Logo\" dialog"), "V-set"));
    await session.step(51, "When the video says \"The V-set profile: its conserved Cys and Trp stand tall\"", () => videoSays(page, "The V-set profile: its conserved Cys and Trp stand tall"));
    await session.step(52, "And user clicks on OK button in \"Profile HMM Logo\" dialog", () => clickOn(page, el("OK button in \"Profile HMM Logo\" dialog")));
    await session.step(53, "Then WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(54, "And \"profile\" property of WebLogo viewer should be \"V-set\"", () => propertyShouldBe(page, "profile", el("WebLogo viewer"), "V-set"));
    await session.step(55, "And the \"positions shown\" reading of WebLogo viewer should be at least 20", () => readingAtLeast(page, "positions shown", el("WebLogo viewer"), 20));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
  });
});
