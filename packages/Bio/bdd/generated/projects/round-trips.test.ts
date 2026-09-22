/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/round-trips.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [GROK-19928, bio.project.round-trip]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {carriesAtLeast} from '../../bindings/annotations.js';
import {knownMonomer, loadedFrom} from '../../bindings/monomer-libs.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnTagLists, columnUnits, columnsEqual} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openDatasetRowsAs, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColors, noBalloons, noErrors, painted, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bio results survive a project save and reopen", () => {
  const session = feature(test, "features/projects/round-trips.feature", import.meta.url);
  test("Sequence Space output survives a save and reopen (GROK-19928)", {tag: ["@serial", "@realizes:GROK-19928", "@realizes:bio.project.round-trip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(26, "Given user opens FASTA_sample dataset keeping the first 20 rows as \"FASTA_sample\"", () => openDatasetRowsAs(page, ds("FASTA_sample"), 20, "FASTA_sample"));
    await session.step(27, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
    await session.step(28, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
    await session.step(29, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(30, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(31, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(32, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(33, "When user saves the current view as project \"bdd-bio-seqspace-{run}\"", () => saveAsProject(page, session.text("bdd-bio-seqspace-{run}")));
    await session.step(34, "And user closes all views", () => closeAllViews(page));
    await session.step(35, "And user opens the \"bdd-bio-seqspace-{run}\" project", () => openProject(page, session.text("bdd-bio-seqspace-{run}")));
    await session.step(36, "Then the table should have 20 rows", () => rowCount(page, 20));
    await session.step(37, "And \"Sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Sequence", "Macromolecule"));
    await session.step(38, "And \"Sequence\" column should have units \"fasta\"", () => columnUnits(page, "Sequence", "fasta"));
    await session.step(39, "And \"Embed_X_1\" column should have no missing values", () => columnComplete(page, "Embed_X_1"));
    await session.step(40, "And \"Embed_Y_1\" column should have no missing values", () => columnComplete(page, "Embed_Y_1"));
    await session.step(41, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(42, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
  test("Antibody numbering survives a save and reopen, and numbers the same again", {tag: ["@serial", "@realizes:GROK-19928", "@realizes:bio.project.round-trip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(47, "Given user opens antibodies dataset keeping the first 20 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 20, "antibodies"));
    await session.step(48, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(49, "And user selects \"kabat\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "kabat", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
    await session.step(50, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(51, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(52, "And a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(53, "When user saves the current view as project \"bdd-bio-numbering-{run}\"", () => saveAsProject(page, session.text("bdd-bio-numbering-{run}")));
    await session.step(54, "And user closes all views", () => closeAllViews(page));
    await session.step(55, "And user opens the \"bdd-bio-numbering-{run}\" project", () => openProject(page, session.text("bdd-bio-numbering-{run}")));
    await session.step(56, "Then the table should have 20 rows", () => rowCount(page, 20));
    await session.step(57, "And \"AntibodyHC\" column should have units \"fasta\"", () => columnUnits(page, "AntibodyHC", "fasta"));
    await session.step(58, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"kabat\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "kabat"));
    await session.step(59, "And the \".positionNames\" tag of \"AntibodyHC (aligned)\" column should list at least 100 values", () => columnTagLists(page, ".positionNames", "AntibodyHC (aligned)", 100));
    await session.step(60, "And \"AntibodyHC (aligned)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned)"));
    await session.step(61, "And \"AntibodyHC\" column should carry at least 7 annotations", () => carriesAtLeast(page, "AntibodyHC", 7));
    await session.step(62, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(63, "And user selects \"kabat\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "kabat", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
    await session.step(64, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(65, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(66, "And a new column \"AntibodyHC (aligned) (2)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned) (2)"));
    await session.step(67, "And \"AntibodyHC (aligned) (2)\" column should hold the same values as \"AntibodyHC (aligned)\" column", () => columnsEqual(page, "AntibodyHC (aligned) (2)", "AntibodyHC (aligned)"));
    await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
  test("A HELM table comes back painted from the same monomer library", {tag: ["@serial", "@realizes:GROK-19928", "@realizes:bio.project.round-trip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(72, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(73, "Then the \"cell type of HELM string\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM string", el("grid"), "helm"));
    await session.step(74, "And the \"cell 1 of HELM string\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of HELM string", el("grid"), 3));
    await session.step(75, "When user saves the current view as project \"bdd-bio-helm-{run}\"", () => saveAsProject(page, session.text("bdd-bio-helm-{run}")));
    await session.step(76, "And user closes all views", () => closeAllViews(page));
    await session.step(77, "And user opens the \"bdd-bio-helm-{run}\" project", () => openProject(page, session.text("bdd-bio-helm-{run}")));
    await session.step(78, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(79, "And \"HELM string\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM string", "Macromolecule"));
    await session.step(80, "And \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
    await session.step(81, "And the \"cell type of HELM string\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM string", el("grid"), "helm"));
    await session.step(82, "And the \"cell 1 of HELM string\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of HELM string", el("grid"), 3));
    await session.step(83, "And the monomer library should be loaded from \"HELMCoreLibrary.json\"", () => loadedFrom(page, "HELMCoreLibrary.json"));
    await session.step(84, "And \"dV\" should be a known \"PEPTIDE\" monomer", () => knownMonomer(page, "dV", "PEPTIDE"));
    await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(86, "And no errors should have been logged", () => noErrors(page));
  });
});
