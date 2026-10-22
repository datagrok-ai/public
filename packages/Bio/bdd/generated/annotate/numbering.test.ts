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
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {carriesAtLeast} from '../../bindings/annotations.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, selectIn, shouldBe, shouldHaveText, shouldHaveValue, shouldOffer, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnTagLists, columnUnits, everyValueMatches, everyValueSameLength, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {removeColumn, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, dialogCloses, noProjectOnServer, openDatasetRowsAs, openProjectWithTable, openSaveDialog, projectsOnServer, saveFromDialog} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Antibody numbering with the bundled immunum engine", () => {
  const session = feature(test, "features/annotate/numbering.feature", import.meta.url);
  test("Antibody numbering with the bundled immunum engine", {tag: ["@journey", "@realizes:bio.annotate.numbering-scheme"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(24, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(25, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await session.step(26, "And \"AntibodyHC\" column should have units \"fasta\"", () => columnUnits(page, "AntibodyHC", "fasta"));
    await run.scenario("The dialog offers the engine and its schemes, IMGT first", async () => {
      await session.step(29, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(30, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(31, "And editor of Sequence input in \"Apply Antibody Numbering\" dialog should have text \"AntibodyHC\"", () => shouldHaveText(page, el("editor of Sequence input in \"Apply Antibody Numbering\" dialog"), "AntibodyHC"));
      await session.step(32, "When user selects \"Immunum\" in Engine input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "Immunum", el("Engine input in \"Apply Antibody Numbering\" dialog")));
      await session.step(33, "Then Engine input in \"Apply Antibody Numbering\" dialog should have value \"Immunum\"", () => shouldHaveValue(page, el("Engine input in \"Apply Antibody Numbering\" dialog"), "Immunum"));
      await session.step(34, "And Scheme input in \"Apply Antibody Numbering\" dialog should offer \"imgt, kabat\"", () => shouldOffer(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt, kabat"));
      await session.step(35, "And Scheme input in \"Apply Antibody Numbering\" dialog should have value \"imgt\"", () => shouldHaveValue(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt"));
    });
    await run.scenario("Kabat numbering aligns the column and annotates its regions", async () => {
      await session.step(38, "When user selects \"kabat\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "kabat", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
      await session.step(39, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
      await session.step(40, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
      await session.step(41, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"kabat\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "kabat"));
      await session.step(42, "And every value of \"AntibodyHC (aligned)\" column should have the same length", () => everyValueSameLength(page, "AntibodyHC (aligned)"));
      await session.step(43, "And \"AntibodyHC (aligned)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned)"));
      await session.step(44, "And every value of \"AntibodyHC (aligned)\" column should match \"[A-Z]{20}\"", () => everyValueMatches(page, "AntibodyHC (aligned)", "[A-Z]{20}"));
      await session.step(45, "And the \".positionNames\" tag of \"AntibodyHC (aligned)\" column should list at least 100 values", () => columnTagLists(page, ".positionNames", "AntibodyHC (aligned)", 100));
      await session.step(46, "And \"AntibodyHC\" column should carry at least 7 annotations", () => carriesAtLeast(page, "AntibodyHC", 7));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The default IMGT numbering aligns the column by IMGT", async () => {
      await session.step(51, "Given user removes \"AntibodyHC (aligned)\" column", () => removeColumn(page, "AntibodyHC (aligned)"));
      await session.step(52, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(53, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(54, "When user selects \"Immunum\" in Engine input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "Immunum", el("Engine input in \"Apply Antibody Numbering\" dialog")));
      await session.step(55, "Then Scheme input in \"Apply Antibody Numbering\" dialog should have value \"imgt\"", () => shouldHaveValue(page, el("Scheme input in \"Apply Antibody Numbering\" dialog"), "imgt"));
      await session.step(56, "When user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
      await session.step(57, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
      await session.step(58, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "imgt"));
      await session.step(59, "And every value of \"AntibodyHC (aligned)\" column should have the same length", () => everyValueSameLength(page, "AntibodyHC (aligned)"));
      await session.step(60, "And \"AntibodyHC (aligned)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned)"));
      await session.step(61, "And the \".positionNames\" tag of \"AntibodyHC (aligned)\" column should list at least 100 values", () => columnTagLists(page, ".positionNames", "AntibodyHC (aligned)", 100));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The numbering comes back from a saved project, and the engine runs again on it", async () => {
      await session.step(66, "Given no project named \"bdd-bio-numbering-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-bio-numbering-{run}")));
      await session.step(67, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(68, "And user types \"bdd-bio-numbering-{run}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("bdd-bio-numbering-{run}"), el("text input in \"Save project\" dialog")));
      await session.step(69, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(70, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(71, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(72, "Then 1 project named \"bdd-bio-numbering-{run}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-bio-numbering-{run}")));
      await session.step(73, "When user closes all views", () => closeAllViews(page));
      await session.step(74, "And user opens the \"bdd-bio-numbering-{run}\" project and waits for its table", () => openProjectWithTable(page, session.text("bdd-bio-numbering-{run}")));
      await session.step(75, "Then the table should have 40 rows", () => rowCount(page, 40));
      await session.step(76, "And \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
      await session.step(77, "And \"AntibodyHC\" column should have units \"fasta\"", () => columnUnits(page, "AntibodyHC", "fasta"));
      await session.step(78, "And the table should have a column \"AntibodyHC (aligned)\"", () => hasColumn(page, "AntibodyHC (aligned)"));
      await session.step(79, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "imgt"));
      await session.step(80, "And every value of \"AntibodyHC (aligned)\" column should have the same length", () => everyValueSameLength(page, "AntibodyHC (aligned)"));
      await session.step(81, "And \"AntibodyHC (aligned)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned)"));
      await session.step(82, "And the \".positionNames\" tag of \"AntibodyHC (aligned)\" column should list at least 100 values", () => columnTagLists(page, ".positionNames", "AntibodyHC (aligned)", 100));
      await session.step(83, "And \"AntibodyHC\" column should carry at least 7 annotations", () => carriesAtLeast(page, "AntibodyHC", 7));
      await session.step(84, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(85, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(86, "When user selects \"Immunum\" in Engine input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "Immunum", el("Engine input in \"Apply Antibody Numbering\" dialog")));
      await session.step(87, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
      await session.step(88, "Then a new column \"AntibodyHC (aligned) (2)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned) (2)"));
      await session.step(89, "And \"AntibodyHC (aligned) (2)\" column should have tag \".numberingScheme\" equal to \"imgt\"", () => columnTag(page, "AntibodyHC (aligned) (2)", ".numberingScheme", "imgt"));
      await session.step(90, "And \"AntibodyHC (aligned) (2)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned) (2)"));
      await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
