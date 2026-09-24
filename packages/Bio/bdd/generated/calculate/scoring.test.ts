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
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnUnits, distinctValues, everyValueBetween, maxInRow, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Identity and similarity scoring", () => {
  const session = feature(test, "features/calculate/scoring.feature", import.meta.url);
  test("Identity and similarity scoring", {tag: ["@journey", "@realizes:bio.calculate.identity", "@realizes:bio.calculate.similarity"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(15, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
    await run.scenario("Identity against the first row", async () => {
      await session.step(18, "When user picks \"Bio > Calculate > Identity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Identity..."));
      await session.step(19, "Then Identity dialog should be visible", () => shouldBe(page, el("Identity dialog"), "visible"));
      await session.step(20, "And editor of Macromolecule input in Identity dialog should have text \"HELM string\"", () => shouldHaveText(page, el("editor of Macromolecule input in Identity dialog"), "HELM string"));
      await session.step(21, "And OK button in Identity dialog should be disabled", () => shouldBe(page, el("OK button in Identity dialog"), "disabled"));
      await session.step(22, "When user enters \"PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0\" into Reference input in Identity dialog", () => enterInto(page, "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0", el("Reference input in Identity dialog")));
      await session.step(23, "Then OK button in Identity dialog should be enabled", () => shouldBe(page, el("OK button in Identity dialog"), "enabled"));
      await session.step(24, "When user clicks on OK button in Identity dialog", () => clickOn(page, el("OK button in Identity dialog")));
      await session.step(25, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(26, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(27, "And a new column \"Identity\" should have been added", () => newColumnNamed(page, "Identity"));
      await session.step(28, "And \"Identity\" column should have no missing values", () => columnComplete(page, "Identity"));
      await session.step(29, "And the value of \"Identity\" column in row 1 should be \"1\"", () => valueInRow(page, "Identity", 1, "1"));
      await session.step(30, "And every value of \"Identity\" column should lie between 0 and 1", () => everyValueBetween(page, "Identity", 0, 1));
      await session.step(31, "And \"Identity\" column should have at least 2 distinct values", () => distinctValues(page, "Identity", 2));
      await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Similarity against the first row peaks on it", async () => {
      await session.step(35, "When user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
      await session.step(36, "Then Similarity dialog should be visible", () => shouldBe(page, el("Similarity dialog"), "visible"));
      await session.step(37, "When user enters \"PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0\" into Reference input in Similarity dialog", () => enterInto(page, "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0", el("Reference input in Similarity dialog")));
      await session.step(38, "And user clicks on OK button in Similarity dialog", () => clickOn(page, el("OK button in Similarity dialog")));
      await session.step(39, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(40, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(41, "And a new column \"Similarity\" should have been added", () => newColumnNamed(page, "Similarity"));
      await session.step(42, "And \"Similarity\" column should have its maximum in row 1", () => maxInRow(page, "Similarity", 1));
      await session.step(43, "And every value of \"Similarity\" column should lie between 0 and 2", () => everyValueBetween(page, "Similarity", 0, 2));
      await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Similarity against the first row scores every row", async () => {
      await session.step(50, "Then \"Similarity\" column should have no missing values", () => columnComplete(page, "Similarity"));
      await session.step(51, "And \"Similarity\" column should have at least 2 distinct values", () => distinctValues(page, "Similarity", 2));
    });
    await run.scenario("Similarity against another reference peaks on that row", async () => {
      await session.step(54, "When user removes \"Similarity\" column", () => removeColumn(page, "Similarity"));
      await session.step(55, "And user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
      await session.step(56, "Then Similarity dialog should be visible", () => shouldBe(page, el("Similarity dialog"), "visible"));
      await session.step(57, "When user enters \"PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$\" into Reference input in Similarity dialog", () => enterInto(page, "PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$", el("Reference input in Similarity dialog")));
      await session.step(58, "And user clicks on OK button in Similarity dialog", () => clickOn(page, el("OK button in Similarity dialog")));
      await session.step(59, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(60, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(61, "And a new column \"Similarity\" should have been added", () => newColumnNamed(page, "Similarity"));
      await session.step(62, "And \"Similarity\" column should have its maximum in row 3", () => maxInRow(page, "Similarity", 3));
    });
    await run.scenario("Similarity leaves no cell blank", async () => {
      await session.step(68, "Then \"Similarity\" column should have no missing values", () => columnComplete(page, "Similarity"));
    });
    run.finish();
  });
});
