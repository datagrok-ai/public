/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/hmm/numbering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [hmm.numbering.engine, hmm.numbering.bio-dialog, hmm.numbering.options, hmm.germlines.top-menu]
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
import {clickOn, followingShouldBe, selectIn, shouldBe, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnTagLists, everyValueSameLength, someValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDatasetRowsAs, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("HMM: ANARCI numbering in Bio's numbering dialog, and germlines from the Bio menu", () => {
  const session = feature(test, "features/hmm/numbering.feature", import.meta.url);
  test("HMM: ANARCI numbering in Bio's numbering dialog, and germlines from the Bio menu", {tag: ["@journey", "@realizes:hmm.numbering.engine", "@realizes:hmm.numbering.bio-dialog", "@realizes:hmm.numbering.options", "@realizes:hmm.germlines.top-menu"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the \"HMM\" package is installed", () => packageInstalled(page, "HMM"));
    await session.step(18, "And user opens antibodies dataset keeping the first 60 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 60, "antibodies"));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(20, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await run.scenario("The dialog offers ANARCI first, with its six schemes and its own options", async () => {
      await session.step(23, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(24, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(25, "And Engine input in \"Apply Antibody Numbering\" dialog should have value \"ANARCI (HMMER)\"", () => shouldHaveValue(page, el("Engine input in \"Apply Antibody Numbering\" dialog"), "ANARCI (HMMER)"));
      await session.step(26, "And Scheme input in \"Apply Antibody Numbering\" dialog should offer \"imgt, kabat, chothia, martin, aho, wolfguy\"", () => shouldOffer(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt, kabat, chothia, martin, aho, wolfguy"));
      await session.step(27, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Species input in \"Apply Antibody Numbering\" dialog"],["\"Bit score threshold\" input in \"Apply Antibody Numbering\" dialog"],["Domain input in \"Apply Antibody Numbering\" dialog"]]), [["Species input in \"Apply Antibody Numbering\" dialog"],["\"Bit score threshold\" input in \"Apply Antibody Numbering\" dialog"],["Domain input in \"Apply Antibody Numbering\" dialog"]]);
    });
    await run.scenario("Chothia numbering of the heavy domain aligns the column", async () => {
      await session.step(33, "When user selects \"chothia\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "chothia", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
      await session.step(34, "And user selects \"heavy\" in Domain input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "heavy", el("Domain input in \"Apply Antibody Numbering\" dialog")));
      await session.step(35, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
      await session.step(36, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
      await session.step(37, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"chothia\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "chothia"));
      await session.step(38, "And every value of \"AntibodyHC (aligned)\" column should have the same length", () => everyValueSameLength(page, "AntibodyHC (aligned)"));
      await session.step(39, "And the \".positionNames\" tag of \"AntibodyHC (aligned)\" column should list at least 100 values", () => columnTagLists(page, ".positionNames", "AntibodyHC (aligned)", 100));
      await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Germlines and Species adds the chain, species and V gene of every chain", async () => {
      await session.step(44, "When user picks \"Bio > Annotate > Germlines and Species (ANARCI)...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Germlines and Species (ANARCI)..."));
      await session.step(45, "Then \"Antibody Germlines (ANARCI)\" dialog should be visible", () => shouldBe(page, el("\"Antibody Germlines (ANARCI)\" dialog"), "visible"));
      await session.step(46, "When user clicks on OK button in \"Antibody Germlines (ANARCI)\" dialog", () => clickOn(page, el("OK button in \"Antibody Germlines (ANARCI)\" dialog")));
      await session.step(47, "Then a new column \"AntibodyHC v gene\" should have been added", () => newColumnNamed(page, "AntibodyHC v gene"));
      await session.step(48, "And \"AntibodyHC v gene\" column should have no missing values", () => columnComplete(page, "AntibodyHC v gene"));
      await session.step(49, "And some value of \"AntibodyHC chain\" column should contain \"Heavy\"", () => someValueContains(page, "AntibodyHC chain", "Heavy"));
      await session.step(50, "And some value of \"AntibodyHC species\" column should contain \"human\"", () => someValueContains(page, "AntibodyHC species", "human"));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
