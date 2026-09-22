/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/similarity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.similarity-search]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, pickFromContextMenu, readingAtLeast, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Similarity Search viewer and its properties", () => {
  const session = feature(test, "features/search/similarity.feature", import.meta.url);
  test("The Similarity Search viewer and its properties", {tag: ["@journey", "@realizes:chem.cp.similarity-search"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("The viewer shows the most similar molecules to the current row", async () => {
      await session.step(16, "When user picks \"Chem > Search > Similarity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Similarity Search..."));
      await session.step(17, "Then Chem Similarity Search viewer should be visible", () => shouldBe(page, el("Chem Similarity Search viewer"), "visible"));
      await session.step(18, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(19, "And the \"metric\" reading of Chem Similarity Search viewer should be \"Tanimoto\"", () => readingReads(page, "metric", el("Chem Similarity Search viewer"), "Tanimoto"));
      await session.step(20, "And the \"fingerprint\" reading of Chem Similarity Search viewer should be \"Morgan\"", () => readingReads(page, "fingerprint", el("Chem Similarity Search viewer"), "Morgan"));
      await session.step(21, "And the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, Morgan\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, Morgan"));
      await session.step(22, "And the \"target row\" reading of Chem Similarity Search viewer should be 1", () => readingIs(page, "target row", el("Chem Similarity Search viewer"), 1));
      await session.step(23, "And the \"scores\" reading of Chem Similarity Search viewer should include the text \"1.00, \"", () => readingIncludes(page, "scores", el("Chem Similarity Search viewer"), "1.00, "));
      await session.step(24, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every fingerprint runs the search again", async () => {
      await session.step(27, "When user picks \"Properties...\" from the context menu of Chem Similarity Search viewer", () => pickFromContextMenu(page, "Properties...", el("Chem Similarity Search viewer")));
      await session.step(28, "And user expands Misc category", () => expand(page, el("Misc category")));
      await session.step(29, "And user selects \"RDKit\" in Fingerprint property", () => selectIn(page, "RDKit", el("Fingerprint property")));
      await session.step(30, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, RDKit\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, RDKit"));
      await session.step(31, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(32, "When user selects \"MACCS\" in Fingerprint property", () => selectIn(page, "MACCS", el("Fingerprint property")));
      await session.step(33, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, MACCS\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, MACCS"));
      await session.step(34, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(35, "When user selects \"AtomPair\" in Fingerprint property", () => selectIn(page, "AtomPair", el("Fingerprint property")));
      await session.step(36, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, AtomPair\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, AtomPair"));
      await session.step(37, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(38, "When user selects \"TopologicalTorsion\" in Fingerprint property", () => selectIn(page, "TopologicalTorsion", el("Fingerprint property")));
      await session.step(39, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, TopologicalTorsion\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, TopologicalTorsion"));
      await session.step(40, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(41, "When user selects \"Morgan\" in Fingerprint property", () => selectIn(page, "Morgan", el("Fingerprint property")));
      await session.step(42, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, Morgan\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, Morgan"));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Limit sets the number of cards", async () => {
      await session.step(47, "When user enters \"5\" into Limit property", () => enterInto(page, "5", el("Limit property")));
      await session.step(48, "Then the \"cards\" reading of Chem Similarity Search viewer should be 5", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 5));
      await session.step(49, "When user enters \"20\" into Limit property", () => enterInto(page, "20", el("Limit property")));
      await session.step(50, "Then the \"cards\" reading of Chem Similarity Search viewer should be 20", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 20));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every metric runs the search again", async () => {
      await session.step(54, "When user selects \"Dice\" in \"Distance Metric\" property", () => selectIn(page, "Dice", el("\"Distance Metric\" property")));
      await session.step(55, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Dice, Morgan\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Dice, Morgan"));
      await session.step(56, "And the \"cards\" reading of Chem Similarity Search viewer should be 20", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 20));
      await session.step(57, "When user selects \"Cosine\" in \"Distance Metric\" property", () => selectIn(page, "Cosine", el("\"Distance Metric\" property")));
      await session.step(58, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Cosine, Morgan\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Cosine, Morgan"));
      await session.step(59, "When user selects \"Tanimoto\" in \"Distance Metric\" property", () => selectIn(page, "Tanimoto", el("\"Distance Metric\" property")));
      await session.step(60, "Then the \"header\" reading of Chem Similarity Search viewer should be \"Tanimoto, Morgan\"", () => readingReads(page, "header", el("Chem Similarity Search viewer"), "Tanimoto, Morgan"));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Size sets the size of every card", async () => {
      await session.step(65, "When user selects \"normal\" in Size property", () => selectIn(page, "normal", el("Size property")));
      await session.step(66, "Then the \"card sizes\" reading of Chem Similarity Search viewer should be \"200x100\"", () => readingReads(page, "card sizes", el("Chem Similarity Search viewer"), "200x100"));
      await session.step(67, "When user selects \"large\" in Size property", () => selectIn(page, "large", el("Size property")));
      await session.step(68, "Then the \"card sizes\" reading of Chem Similarity Search viewer should be \"300x150\"", () => readingReads(page, "card sizes", el("Chem Similarity Search viewer"), "300x150"));
      await session.step(69, "When user selects \"small\" in Size property", () => selectIn(page, "small", el("Size property")));
      await session.step(70, "Then the \"card sizes\" reading of Chem Similarity Search viewer should be \"120x60\"", () => readingReads(page, "card sizes", el("Chem Similarity Search viewer"), "120x60"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Molecule Properties adds the chosen columns to the cards", async () => {
      await session.step(74, "When user clicks on \"...\" button in \"Molecule Properties\" property", () => clickOn(page, el("\"...\" button in \"Molecule Properties\" property")));
      await session.step(75, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(76, "When user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(77, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(78, "Then the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"molregno\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "molregno"));
      await session.step(79, "And the \"text of cell 3 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"NumSaturatedHeterocycles\"", () => readingReads(page, "text of cell 3 of __name", el("grid viewer in \"Select columns...\" dialog"), "NumSaturatedHeterocycles"));
      await session.step(80, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(81, "Then the \"card properties\" reading of Chem Similarity Search viewer should be \"molregno, NumSaturatedHeterocycles\"", () => readingReads(page, "card properties", el("Chem Similarity Search viewer"), "molregno, NumSaturatedHeterocycles"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cutoff 1 leaves the molecules that score 1", async () => {
      await session.step(85, "When user enters \"1\" into Cutoff property", () => enterInto(page, "1", el("Cutoff property")));
      await session.step(86, "Then the \"min score\" reading of Chem Similarity Search viewer should be 1", () => readingIs(page, "min score", el("Chem Similarity Search viewer"), 1));
      await session.step(87, "And the \"cards\" reading of Chem Similarity Search viewer should be at least 1", () => readingAtLeast(page, "cards", el("Chem Similarity Search viewer"), 1));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
