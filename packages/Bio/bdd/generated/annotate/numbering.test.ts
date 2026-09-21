/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotate/numbering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.annotate.numbering-scheme]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {carriesAtLeast} from '../../bindings/annotations.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnUnits, everyValueSameLength} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {callWith, resultTableColumns, resultTableFilled} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDatasetRowsAs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Antibody numbering with the bundled immunum engine", () => {
  const session = feature(test, "features/annotate/numbering.feature", import.meta.url);
  test("Antibody numbering with the bundled immunum engine", {tag: ["@journey", "@realizes:bio.annotate.numbering-scheme"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(11, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(12, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await session.step(13, "And \"AntibodyHC\" column should have units \"fasta\"", () => columnUnits(page, "AntibodyHC", "fasta"));
    await run.scenario("The dialog offers the engine and its schemes, IMGT first", async () => {
      await session.step(16, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(17, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(18, "And editor of Sequence input in \"Apply Antibody Numbering\" dialog should have text \"AntibodyHC\"", () => shouldHaveText(page, el("editor of Sequence input in \"Apply Antibody Numbering\" dialog"), "AntibodyHC"));
      await session.step(19, "And Scheme input in \"Apply Antibody Numbering\" dialog should have value \"imgt\"", () => shouldHaveValue(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt"));
    });
    await run.scenario("Kabat numbering aligns the column and annotates its regions", async () => {
      await session.step(22, "When user selects \"kabat\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "kabat", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
      await session.step(23, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
      await session.step(24, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(25, "And a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
      await session.step(26, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"kabat\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "kabat"));
      await session.step(27, "And every value of \"AntibodyHC (aligned)\" column should have the same length", () => everyValueSameLength(page, "AntibodyHC (aligned)"));
      await session.step(28, "And \"AntibodyHC (aligned)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned)"));
      await session.step(29, "And \"AntibodyHC\" column should carry at least 7 annotations", () => carriesAtLeast(page, "AntibodyHC", 7));
      await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The engine's own result honours the five-column contract", async () => {
      await session.step(34, "When user calls \"Bio:immunumAntibodyNumbering\" function with:", () => callWith(page, "Bio:immunumAntibodyNumbering", [["df","table"],["seqCol","column:AntibodyHC"],["scheme","imgt"]]));
      await session.step(38, "Then the result should be a table with columns \"position_names, chain_type, annotations_json, numbering_detail, numbering_map\"", () => resultTableColumns(page, "position_names, chain_type, annotations_json, numbering_detail, numbering_map"));
      await session.step(39, "And every column of the result table should be filled in row 1", () => resultTableFilled(page, 1));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
