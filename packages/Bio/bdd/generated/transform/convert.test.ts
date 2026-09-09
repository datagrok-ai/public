/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/convert.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.transform.convert-notation, bio.calculate.extract-region, bio.transform.split-to-monomers]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, columnUnits, everyValueMatches, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Transforming a fasta column", () => {
  const session = feature(test, "features/transform/convert.feature", import.meta.url);
  test("Transforming a fasta column", {tag: ["@journey", "@realizes:bio.transform.convert-notation", "@realizes:bio.calculate.extract-region", "@realizes:bio.transform.split-to-monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens filter_FASTA dataset keeping the first 9 rows", () => openDatasetRows(page, ds("filter_FASTA"), 9));
    await session.step(12, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(13, "Then the table should have 9 rows", () => rowCount(page, 9));
    await session.step(14, "And \"fasta\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "fasta", "Macromolecule"));
    await session.step(15, "And \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await session.step(16, "And \"fasta\" column should have tag \"alphabet\" equal to \"PT\"", () => columnTag(page, "fasta", "alphabet", "PT"));
    await run.scenario("Extract Region proposes the whole range and cuts it out as a sequence column", async () => {
      await session.step(19, "When user picks \"Bio > Calculate > Extract Region...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Extract Region..."));
      await session.step(20, "Then \"Get Sequence Region\" dialog should be visible", () => shouldBe(page, el("\"Get Sequence Region\" dialog"), "visible"));
      await session.step(21, "And Start input in \"Get Sequence Region\" dialog should have value \"1\"", () => shouldHaveValue(page, el("Start input in \"Get Sequence Region\" dialog"), "1"));
      await session.step(22, "And End input in \"Get Sequence Region\" dialog should have value \"38\"", () => shouldHaveValue(page, el("End input in \"Get Sequence Region\" dialog"), "38"));
      await session.step(23, "And \"Column name\" input in \"Get Sequence Region\" dialog should have value \"fasta: (1-38)\"", () => shouldHaveValue(page, el("\"Column name\" input in \"Get Sequence Region\" dialog"), "fasta: (1-38)"));
      await session.step(24, "When user selects \"3\" in Start input in \"Get Sequence Region\" dialog", () => selectIn(page, "3", el("Start input in \"Get Sequence Region\" dialog")));
      await session.step(25, "And user selects \"6\" in End input in \"Get Sequence Region\" dialog", () => selectIn(page, "6", el("End input in \"Get Sequence Region\" dialog")));
      await session.step(26, "And user enters \"region 3-6\" into \"Column name\" input in \"Get Sequence Region\" dialog", () => enterInto(page, "region 3-6", el("\"Column name\" input in \"Get Sequence Region\" dialog")));
      await session.step(27, "And user clicks on OK button in \"Get Sequence Region\" dialog", () => clickOn(page, el("OK button in \"Get Sequence Region\" dialog")));
      await session.step(28, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(29, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(30, "And a new column \"region 3-6\" should have been added", () => newColumnNamed(page, "region 3-6"));
      await session.step(31, "And \"region 3-6\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "region 3-6", "Macromolecule"));
      await session.step(32, "And \"region 3-6\" column should have units \"fasta\"", () => columnUnits(page, "region 3-6", "fasta"));
      await session.step(33, "And every value of \"region 3-6\" column should match \"^[A-Z]{4}$\"", () => everyValueMatches(page, "region 3-6", "^[A-Z]{4}$"));
      await session.step(34, "And the value of \"region 3-6\" column in row 1 should be \"YKET\"", () => valueInRow(page, "region 3-6", 1, "YKET"));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Convert Sequence Notation proposes separator and writes the column in it", async () => {
      await session.step(38, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
      await session.step(39, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
      await session.step(40, "And \"Convert Sequence Notation\" dialog should contain text \"Current notation: fasta\"", () => shouldContainText(page, el("\"Convert Sequence Notation\" dialog"), "Current notation: fasta"));
      await session.step(41, "And \"Convert to\" input in \"Convert Sequence Notation\" dialog should have value \"separator\"", () => shouldHaveValue(page, el("\"Convert to\" input in \"Convert Sequence Notation\" dialog"), "separator"));
      await session.step(42, "And Separator input in \"Convert Sequence Notation\" dialog should have value \"-\"", () => shouldHaveValue(page, el("Separator input in \"Convert Sequence Notation\" dialog"), "-"));
      await session.step(43, "When user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
      await session.step(44, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(45, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(46, "And a new column matching \"^separator\\(fasta\\)\" should have been added", () => newColumnMatching(page, "^separator\\(fasta\\)"));
      await session.step(47, "And \"separator(fasta)\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "separator(fasta)", "Macromolecule"));
      await session.step(48, "And \"separator(fasta)\" column should have units \"separator\"", () => columnUnits(page, "separator(fasta)", "separator"));
      await session.step(49, "And \"separator(fasta)\" column should have tag \"separator\" equal to \"-\"", () => columnTag(page, "separator(fasta)", "separator", "-"));
      await session.step(50, "And the value of \"separator(fasta)\" column in row 1 should be \"M-D-Y-K-E-T-L-L-M-P-K-T-D-F-P-M-R-G-G-L-P-N-K-E-P-Q-I-Q-E-K-W\"", () => valueInRow(page, "separator(fasta)", 1, "M-D-Y-K-E-T-L-L-M-P-K-T-D-F-P-M-R-G-G-L-P-N-K-E-P-Q-I-Q-E-K-W"));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Converting to HELM wraps every sequence in a PEPTIDE polymer", async () => {
      await session.step(54, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
      await session.step(55, "And user selects \"helm\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "helm", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
      await session.step(56, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
      await session.step(57, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(58, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(59, "And a new column \"helm(fasta)\" should have been added", () => newColumnNamed(page, "helm(fasta)"));
      await session.step(60, "And \"helm(fasta)\" column should have units \"helm\"", () => columnUnits(page, "helm(fasta)", "helm"));
      await session.step(61, "And every value of \"helm(fasta)\" column should match \"^PEPTIDE1\\{([A-Z]\\.)*[A-Z]\\}\\$\"", () => everyValueMatches(page, "helm(fasta)", "^PEPTIDE1\\{([A-Z]\\.)*[A-Z]\\}\\$"));
      await session.step(62, "And the value of \"helm(fasta)\" column in row 1 should be \"PEPTIDE1{M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W}$$$$\"", () => valueInRow(page, "helm(fasta)", 1, "PEPTIDE1{M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W}$$$$"));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Split to Monomers adds one Monomer column per position of the longest sequence", async () => {
      await session.step(66, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
      await session.step(67, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
      await session.step(68, "And editor of Sequence input in \"Split to Monomers\" dialog should have text \"fasta\"", () => shouldHaveText(page, el("editor of Sequence input in \"Split to Monomers\" dialog"), "fasta"));
      await session.step(69, "When user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
      await session.step(70, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(71, "And 38 new columns should have been added", () => newColumnsCount(page, 38));
      await session.step(72, "And a new column \"1\" should have been added", () => newColumnNamed(page, "1"));
      await session.step(73, "And a new column \"38\" should have been added", () => newColumnNamed(page, "38"));
      await session.step(74, "And \"1\" column should have semantic type \"Monomer\"", () => columnSemType(page, "1", "Monomer"));
      await session.step(75, "And the value of \"1\" column in row 1 should be \"M\"", () => valueInRow(page, "1", 1, "M"));
      await session.step(76, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
