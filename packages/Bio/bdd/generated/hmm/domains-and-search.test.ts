/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/hmm/domains-and-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [hmm.domains.top-menu, hmm.search.top-menu, hmm.panels.anarci, hmm.panels.pfam]
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
import {clickOn, expand, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, someValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelOpen, openDatasetRowsAs, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("HMM: protein domains, profile search and the sequence panels", () => {
  const session = feature(test, "features/hmm/domains-and-search.feature", import.meta.url);
  test("HMM: protein domains, profile search and the sequence panels", {tag: ["@journey", "@realizes:hmm.domains.top-menu", "@realizes:hmm.search.top-menu", "@realizes:hmm.panels.anarci", "@realizes:hmm.panels.pfam"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(20, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await run.scenario("An antibody chain shows its ANARCI and Pfam Domains panes", async () => {
      await session.step(25, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(26, "When user clicks on the \"cell 2 of AntibodyHC\" area of grid", () => clickArea(page, "cell 2 of AntibodyHC", el("grid")));
      await session.step(27, "Then the \"current row\" reading of grid should be 2", () => readingIs(page, "current row", el("grid"), 2));
      await session.step(28, "When user expands Bioinformatics accordion header in context panel", () => expand(page, el("Bioinformatics accordion header in context panel")));
      await session.step(29, "And user expands ANARCI accordion header in context panel", () => expand(page, el("ANARCI accordion header in context panel")));
      await session.step(30, "Then \"ANARCI\" pane in context panel should contain the text \"Heavy\"", () => shouldContainText(page, el("\"ANARCI\" pane in context panel"), "Heavy"));
      await session.step(31, "And \"ANARCI\" pane in context panel should contain the text \"IGHV\"", () => shouldContainText(page, el("\"ANARCI\" pane in context panel"), "IGHV"));
      await session.step(32, "And \"ANARCI\" pane in context panel should contain the text \"CDR3 (IMGT)\"", () => shouldContainText(page, el("\"ANARCI\" pane in context panel"), "CDR3 (IMGT)"));
      await session.step(33, "When user expands \"Pfam Domains\" accordion header in context panel", () => expand(page, el("\"Pfam Domains\" accordion header in context panel")));
      await session.step(34, "Then \"Pfam Domains\" pane in context panel should contain the text \"V-set\"", () => shouldContainText(page, el("\"Pfam Domains\" pane in context panel"), "V-set"));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Find Domains names the immunoglobulin domains of the chains", async () => {
      await session.step(39, "When user picks \"Bio > Annotate > Find Domains (HMMER)...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Find Domains (HMMER)..."));
      await session.step(40, "Then \"Find Domains (HMMER)\" dialog should be visible", () => shouldBe(page, el("\"Find Domains (HMMER)\" dialog"), "visible"));
      await session.step(41, "And editor of Sequence input in \"Find Domains (HMMER)\" dialog should have text \"AntibodyHC\"", () => shouldHaveText(page, el("editor of Sequence input in \"Find Domains (HMMER)\" dialog"), "AntibodyHC"));
      await session.step(42, "When user clicks on OK button in \"Find Domains (HMMER)\" dialog", () => clickOn(page, el("OK button in \"Find Domains (HMMER)\" dialog")));
      await session.step(43, "Then a new column \"AntibodyHC domains\" should have been added", () => newColumnNamed(page, "AntibodyHC domains"));
      await session.step(44, "And some value of \"AntibodyHC domains\" column should contain \"V-set\"", () => someValueContains(page, "AntibodyHC domains", "V-set"));
      await session.step(45, "And some value of \"AntibodyHC domains\" column should contain \"C1-set\"", () => someValueContains(page, "AntibodyHC domains", "C1-set"));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Profile HMM Search with the built-in V-set family adds its columns", async () => {
      await session.step(50, "When user picks \"Bio > Search > Profile HMM Search...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Profile HMM Search..."));
      await session.step(51, "Then \"Profile HMM Search\" dialog should be visible", () => shouldBe(page, el("\"Profile HMM Search\" dialog"), "visible"));
      await session.step(52, "When user clicks on OK button in \"Profile HMM Search\" dialog", () => clickOn(page, el("OK button in \"Profile HMM Search\" dialog")));
      await session.step(53, "Then a new column \"V-set significant\" should have been added", () => newColumnNamed(page, "V-set significant"));
      await session.step(54, "And a new column \"V-set E-value\" should have been added", () => newColumnNamed(page, "V-set E-value"));
      await session.step(55, "And a new column \"V-set score\" should have been added", () => newColumnNamed(page, "V-set score"));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
