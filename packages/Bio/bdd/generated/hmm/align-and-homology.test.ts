/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/hmm/align-and-homology.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [hmm.align.msa-engine, hmm.align.bio-dialog, hmm.homology.top-menu, hmm.logo.weblogo, hmm.build.top-menu]
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
import {clickOn, enterInto, selectIn, shouldBe, shouldHaveValue, shouldOfferChoice} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, columnUnits, everyValueSameLength, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectFirstRows, selectNoRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDatasetRowsAs, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingAtLeast} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("HMM: profile alignment, homology search, profile logos and built profiles from the Bio menu", () => {
  const session = feature(test, "features/hmm/align-and-homology.feature", import.meta.url);
  test("HMM: profile alignment, homology search, profile logos and built profiles from the Bio menu", {tag: ["@journey", "@realizes:hmm.align.msa-engine", "@realizes:hmm.align.bio-dialog", "@realizes:hmm.homology.top-menu", "@realizes:hmm.logo.weblogo", "@realizes:hmm.build.top-menu"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(21, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(22, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(23, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await run.scenario("The MSA dialog offers hmmalign beside Kalign and aligns the chains to the germline profile", async () => {
      await session.step(26, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
      await session.step(27, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
      await session.step(28, "And Engine input in MSA dialog should have value \"Kalign\"", () => shouldHaveValue(page, el("Engine input in MSA dialog"), "Kalign"));
      await session.step(29, "And Engine input in MSA dialog should offer the choice \"Profile HMM (hmmalign)\"", () => shouldOfferChoice(page, el("Engine input in MSA dialog"), "Profile HMM (hmmalign)"));
      await session.step(30, "When user selects \"Profile HMM (hmmalign)\" in Engine input in MSA dialog", () => selectIn(page, "Profile HMM (hmmalign)", el("Engine input in MSA dialog")));
      await session.step(31, "And user selects \"Antibody/TCR germline (auto)\" in Profile input in MSA dialog", () => selectIn(page, "Antibody/TCR germline (auto)", el("Profile input in MSA dialog")));
      await session.step(32, "And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
      await session.step(33, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
      await session.step(34, "And a new column \"msa(AntibodyHC)\" should have been added", () => newColumnNamed(page, "msa(AntibodyHC)"));
      await session.step(35, "And \"msa(AntibodyHC)\" column should have units \"fasta\"", () => columnUnits(page, "msa(AntibodyHC)", "fasta"));
      await session.step(36, "And \"msa(AntibodyHC)\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "msa(AntibodyHC)", "aligned", "SEQ.MSA"));
      await session.step(37, "And every value of \"msa(AntibodyHC)\" column should have the same length", () => everyValueSameLength(page, "msa(AntibodyHC)"));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Homology Search takes the current row as its query", async () => {
      await session.step(42, "When user makes row 3 current", () => makeRowCurrent(page, 3));
      await session.step(43, "And user picks \"Bio > Search > Homology Search (phmmer)...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Homology Search (phmmer)..."));
      await session.step(44, "Then \"Homology Search (phmmer)\" dialog should be visible", () => shouldBe(page, el("\"Homology Search (phmmer)\" dialog"), "visible"));
      await session.step(45, "And Row input in \"Homology Search (phmmer)\" dialog should have value \"3\"", () => shouldHaveValue(page, el("Row input in \"Homology Search (phmmer)\" dialog"), "3"));
      await session.step(46, "When user clicks on OK button in \"Homology Search (phmmer)\" dialog", () => clickOn(page, el("OK button in \"Homology Search (phmmer)\" dialog")));
      await session.step(47, "Then a new column \"phmmer row 3 E-value\" should have been added", () => newColumnNamed(page, "phmmer row 3 E-value"));
      await session.step(48, "And a new column \"phmmer row 3 score\" should have been added", () => newColumnNamed(page, "phmmer row 3 score"));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Profile HMM Logo docks a WebLogo of a family of the built-in library", async () => {
      await session.step(53, "When user picks \"Bio > Analyze > Profile HMM Logo...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Profile HMM Logo..."));
      await session.step(54, "Then \"Profile HMM Logo\" dialog should be visible", () => shouldBe(page, el("\"Profile HMM Logo\" dialog"), "visible"));
      await session.step(55, "When user selects \"Pfam family (built-in)\" in \"Profile from\" input in \"Profile HMM Logo\" dialog", () => selectIn(page, "Pfam family (built-in)", el("\"Profile from\" input in \"Profile HMM Logo\" dialog")));
      await session.step(56, "And user selects \"C1-set\" in Family input in \"Profile HMM Logo\" dialog", () => selectIn(page, "C1-set", el("Family input in \"Profile HMM Logo\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Profile HMM Logo\" dialog", () => clickOn(page, el("OK button in \"Profile HMM Logo\" dialog")));
      await session.step(58, "Then WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
      await session.step(59, "And the \"positions shown\" reading of WebLogo viewer should be at least 20", () => readingAtLeast(page, "positions shown", el("WebLogo viewer"), 20));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Build Profile HMM from the selected rows searches the table and shows a logo whose letters select rows", async () => {
      await session.step(64, "When user selects the first 15 rows", () => selectFirstRows(page, 15));
      await session.step(65, "And user picks \"Bio > Search > Build Profile HMM...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Build Profile HMM..."));
      await session.step(66, "Then \"Build Profile HMM\" dialog should be visible", () => shouldBe(page, el("\"Build Profile HMM\" dialog"), "visible"));
      await session.step(67, "And Rows input in \"Build Profile HMM\" dialog should have value \"Selected rows\"", () => shouldHaveValue(page, el("Rows input in \"Build Profile HMM\" dialog"), "Selected rows"));
      await session.step(68, "When user enters \"heavy15\" into Name input in \"Build Profile HMM\" dialog", () => enterInto(page, "heavy15", el("Name input in \"Build Profile HMM\" dialog")));
      await session.step(69, "And user clicks on OK button in \"Build Profile HMM\" dialog", () => clickOn(page, el("OK button in \"Build Profile HMM\" dialog")));
      await session.step(70, "Then a new column \"heavy15 E-value\" should have been added", () => newColumnNamed(page, "heavy15 E-value"));
      await session.step(72, "And second WebLogo viewer should be visible", () => shouldBe(page, el("second WebLogo viewer"), "visible"));
      await session.step(73, "When user selects no rows", () => selectNoRows(page));
      await session.step(75, "And user clicks on the \"monomer V at position 2\" area of second WebLogo viewer", () => clickArea(page, "monomer V at position 2", el("second WebLogo viewer")), undefined, "Click the valine at position 2 of the profile");
      await session.step(76, "Then 35 rows should be selected", () => selectedRowCount(page, 35));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
