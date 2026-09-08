/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/scoring.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.calculate.identity, bio.calculate.similarity]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {alignmentLength, bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnIncomplete, columnUnits, distinctValues, everyValueBetween, hasColumn, maxInRow, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterNotNull, filterPasses, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {callWith, resultEmpty} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Identity and similarity scoring", () => {
  const session = feature(test, "features/calculate/scoring.feature", import.meta.url);
  test("Identity and similarity scoring", {tag: ["@journey", "@realizes:bio.calculate.identity", "@realizes:bio.calculate.similarity"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
    await run.scenario("Identity against the first row", async () => {
      await session.step(14, "When user picks \"Bio > Calculate > Identity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Identity..."));
      await session.step(15, "Then Identity dialog should be visible", () => shouldBe(page, el("Identity dialog"), "visible"));
      await session.step(16, "And editor of Macromolecule input in Identity dialog should have text \"HELM string\"", () => shouldHaveText(page, el("editor of Macromolecule input in Identity dialog"), "HELM string"));
      await session.step(17, "And OK button in Identity dialog should be disabled", () => shouldBe(page, el("OK button in Identity dialog"), "disabled"));
      await session.step(18, "When user enters \"PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0\" into Reference input in Identity dialog", () => enterInto(page, "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0", el("Reference input in Identity dialog")));
      await session.step(19, "Then OK button in Identity dialog should be enabled", () => shouldBe(page, el("OK button in Identity dialog"), "enabled"));
      await session.step(20, "When user clicks on OK button in Identity dialog", () => clickOn(page, el("OK button in Identity dialog")));
      await session.step(21, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(22, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(23, "And a new column \"Identity\" should have been added", () => newColumnNamed(page, "Identity"));
      await session.step(24, "And \"Identity\" column should have no missing values", () => columnComplete(page, "Identity"));
      await session.step(25, "And the value of \"Identity\" column in row 1 should be \"1\"", () => valueInRow(page, "Identity", 1, "1"));
      await session.step(26, "And every value of \"Identity\" column should lie between 0 and 1", () => everyValueBetween(page, "Identity", 0, 1));
      await session.step(27, "And \"Identity\" column should have at least 2 distinct values", () => distinctValues(page, "Identity", 2));
      await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Similarity against the first row peaks on it", async () => {
      await session.step(31, "When user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
      await session.step(32, "Then Similarity dialog should be visible", () => shouldBe(page, el("Similarity dialog"), "visible"));
      await session.step(33, "When user enters \"PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0\" into Reference input in Similarity dialog", () => enterInto(page, "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0", el("Reference input in Similarity dialog")));
      await session.step(34, "And user clicks on OK button in Similarity dialog", () => clickOn(page, el("OK button in Similarity dialog")));
      await session.step(35, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(36, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(37, "And a new column \"Similarity\" should have been added", () => newColumnNamed(page, "Similarity"));
      await session.step(38, "And \"Similarity\" column should have its maximum in row 1", () => maxInRow(page, "Similarity", 1));
      await session.step(39, "And every value of \"Similarity\" column should lie between 0 and 2", () => everyValueBetween(page, "Similarity", 0, 2));
      await session.step(40, "And the table should have a column \"Identity\"", () => hasColumn(page, "Identity"));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And \"Similarity\" column should have missing values", () => columnIncomplete(page, "Similarity"));
      await session.step(44, "When user removes \"Similarity\" column", () => removeColumn(page, "Similarity"));
      await session.step(45, "And user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
      await session.step(46, "Then Similarity dialog should be visible", () => shouldBe(page, el("Similarity dialog"), "visible"));
      await session.step(47, "When user enters \"PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$\" into Reference input in Similarity dialog", () => enterInto(page, "PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$", el("Reference input in Similarity dialog")));
      await session.step(48, "And user clicks on OK button in Similarity dialog", () => clickOn(page, el("OK button in Similarity dialog")));
      await session.step(49, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(50, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(51, "And a new column \"Similarity\" should have been added", () => newColumnNamed(page, "Similarity"));
      await session.step(52, "And \"Similarity\" column should have its maximum in row 3", () => maxInRow(page, "Similarity", 3));
      await session.step(53, "When user filters rows where \"Similarity\" is not null", () => filterNotNull(page, "Similarity"));
      await session.step(54, "Then 2 row should pass the filter", () => filterPasses(page, 2));
    });
    await run.scenario("The scoring functions answer an empty sequence with nothing, not an error", async () => {
      await session.step(57, "When user calls \"Bio:seqIdentity\" function with:", () => callWith(page, "Bio:seqIdentity", [["seq",""],["ref","PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0"]]));
      await session.step(60, "Then the result should be empty", () => resultEmpty(page));
      await session.step(61, "When user calls \"Bio:sequenceAlignment\" function with:", () => callWith(page, "Bio:sequenceAlignment", [["alignType","Global alignment"],["alignTable","BLOSUM62"],["gap","-10"],["seq1","MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"],["seq2","MIEVFLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL"]]));
      await session.step(67, "Then the result should be an alignment of at least 37 positions", () => alignmentLength(page, 37));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
