/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/msa.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.msa]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, selectIn, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnType, columnUnits, everyValueMatches, everyValueSameLength, sameLengthPerGroup} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCalculated} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Multiple sequence alignment with kalign", () => {
  const session = feature(test, "features/analyze/msa.feature", import.meta.url);
  test("Multiple sequence alignment with kalign", {tag: ["@journey", "@realizes:bio.analyze.msa"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens filter_FASTA dataset keeping the first 9 rows", () => openDatasetRows(page, ds("filter_FASTA"), 9));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await session.step(12, "And \"fasta\" column should have tag \"alphabet\" equal to \"PT\"", () => columnTag(page, "fasta", "alphabet", "PT"));
    await run.scenario("The dialog opens in kalign mode on the sequence column", async () => {
      await session.step(15, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
      await session.step(16, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
      await session.step(17, "And editor of Sequence input in MSA dialog should have text \"fasta\"", () => shouldHaveText(page, el("editor of Sequence input in MSA dialog"), "fasta"));
      await session.step(18, "And Clusters input in MSA dialog should be visible", () => shouldBe(page, el("Clusters input in MSA dialog"), "visible"));
      await session.step(19, "And Engine input in MSA dialog should be hidden", () => shouldBe(page, el("Engine input in MSA dialog"), "hidden"));
      await session.step(20, "And MSA dialog should contain text \"Kalign version\"", () => shouldContainText(page, el("MSA dialog"), "Kalign version"));
      await session.step(21, "And \"Selected Rows Only\" checkbox in MSA dialog should be unchecked", () => shouldBe(page, el("\"Selected Rows Only\" checkbox in MSA dialog"), "unchecked"));
    });
    await run.scenario("Alignment parameters toggles the kalign penalties", async () => {
      await session.step(24, "Then \"Gap open\" input in MSA dialog should be hidden", () => shouldBe(page, el("\"Gap open\" input in MSA dialog"), "hidden"));
      await session.step(25, "When user clicks on \"Alignment parameters\" button in MSA dialog", () => clickOn(page, el("\"Alignment parameters\" button in MSA dialog")));
      await session.step(26, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Gap open\" input in MSA dialog"],["\"Gap extend\" input in MSA dialog"],["\"Terminal gap\" input in MSA dialog"]]));
      await session.step(30, "When user clicks on \"Alignment parameters\" button in MSA dialog", () => clickOn(page, el("\"Alignment parameters\" button in MSA dialog")));
      await session.step(31, "Then the following elements should be hidden:", () => followingShouldBe(page, "hidden", [["\"Gap open\" input in MSA dialog"],["\"Gap extend\" input in MSA dialog"],["\"Terminal gap\" input in MSA dialog"]]));
    });
    await run.scenario("OK aligns the column", async () => {
      await session.step(37, "When user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
      await session.step(38, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
      await session.step(39, "And a new column \"msa(fasta)\" should have been added", () => newColumnNamed(page, "msa(fasta)"));
      await session.step(40, "And \"msa(fasta)\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "msa(fasta)", "Macromolecule"));
      await session.step(41, "And \"msa(fasta)\" column should have units \"fasta\"", () => columnUnits(page, "msa(fasta)", "fasta"));
      await session.step(42, "And \"msa(fasta)\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "msa(fasta)", "aligned", "SEQ.MSA"));
      await session.step(43, "And \"msa(fasta)\" column should have no missing values", () => columnComplete(page, "msa(fasta)"));
      await session.step(44, "And every value of \"msa(fasta)\" column should have the same length", () => everyValueSameLength(page, "msa(fasta)"));
      await session.step(45, "And every value of \"msa(fasta)\" column should match \"^[A-Z-]+$\"", () => everyValueMatches(page, "msa(fasta)", "^[A-Z-]+$"));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cluster column aligns each cluster on its own", async () => {
      await session.step(50, "When user adds a calculated column \"Clusters\" with formula \"Length(${fasta}) % 2\"", () => addCalculated(page, "Clusters", "Length(${fasta}) % 2"));
      await session.step(51, "Then \"Clusters\" column should have type \"int\"", () => columnType(page, "Clusters", "int"));
      await session.step(52, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
      await session.step(53, "And user selects \"Clusters\" in Clusters input in MSA dialog", () => selectIn(page, "Clusters", el("Clusters input in MSA dialog")));
      await session.step(54, "And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
      await session.step(55, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
      await session.step(56, "And a new column \"msa(fasta) (2)\" should have been added", () => newColumnNamed(page, "msa(fasta) (2)"));
      await session.step(57, "And \"msa(fasta) (2)\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "msa(fasta) (2)", "aligned", "SEQ.MSA"));
      await session.step(58, "And \"msa(fasta) (2)\" column should have no missing values", () => columnComplete(page, "msa(fasta) (2)"));
      await session.step(59, "And every value of \"msa(fasta) (2)\" column should have the same length within each \"Clusters\" value", () => sameLengthPerGroup(page, "msa(fasta) (2)", "Clusters"));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
