/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/round-trips.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.project.round-trip]
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
import {clickOn, selectIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnTag, columnTagLists, columnUnits, columnsEqual} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDatasetRowsAs, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bio results survive a project save and reopen", () => {
  const session = feature(test, "features/projects/round-trips.feature", import.meta.url);
  test("Antibody numbering survives a save and reopen, and numbers the same again", {tag: ["@realizes:bio.project.round-trip"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(27, "Given user opens antibodies dataset keeping the first 20 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 20, "antibodies"));
    await session.step(28, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(29, "And user selects \"kabat\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "kabat", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
    await session.step(30, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(31, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
    await session.step(32, "When user saves the current view as project \"bdd-bio-numbering-{run}\"", () => saveAsProject(page, session.text("bdd-bio-numbering-{run}")));
    await session.step(33, "And user closes all views", () => closeAllViews(page));
    await session.step(34, "And user opens the \"bdd-bio-numbering-{run}\" project", () => openProject(page, session.text("bdd-bio-numbering-{run}")));
    await session.step(35, "Then the table should have 20 rows", () => rowCount(page, 20));
    await session.step(36, "And \"AntibodyHC\" column should have units \"fasta\"", () => columnUnits(page, "AntibodyHC", "fasta"));
    await session.step(37, "And \"AntibodyHC (aligned)\" column should have tag \".numberingScheme\" equal to \"kabat\"", () => columnTag(page, "AntibodyHC (aligned)", ".numberingScheme", "kabat"));
    await session.step(38, "And the \".positionNames\" tag of \"AntibodyHC (aligned)\" column should list at least 100 values", () => columnTagLists(page, ".positionNames", "AntibodyHC (aligned)", 100));
    await session.step(39, "And \"AntibodyHC (aligned)\" column should have no missing values", () => columnComplete(page, "AntibodyHC (aligned)"));
    await session.step(40, "And \"AntibodyHC\" column should carry at least 7 annotations", () => carriesAtLeast(page, "AntibodyHC", 7));
    await session.step(41, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
    await session.step(42, "And user selects \"kabat\" in Scheme input in \"Apply Antibody Numbering\" dialog", () => selectIn(page, "kabat", el("Scheme input in \"Apply Antibody Numbering\" dialog")));
    await session.step(43, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
    await session.step(44, "Then a new column \"AntibodyHC (aligned) (2)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned) (2)"));
    await session.step(45, "And \"AntibodyHC (aligned) (2)\" column should hold the same values as \"AntibodyHC (aligned)\" column", () => columnsEqual(page, "AntibodyHC (aligned) (2)", "AntibodyHC (aligned)"));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
});
