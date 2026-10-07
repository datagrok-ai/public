/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/crux-bio-hosts.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {rgroupLabels} from '../../bindings/monomer-form.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {call} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {autostartsCompleted, openDataset, packageInstalled, sketcherIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, doubleClickArea, hasArea, pickFromOpenMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Crux in Bio's hosts", () => {
  const session = feature(test, "features/sketcher/crux-bio-hosts.feature", import.meta.url);
  test("A 339-atom peptide of the atomic-level demo opens whole in Crux, and an edit there lands", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "When user calls \"Bio:demoBioAtomicLevel\" function", () => call(page, "Bio:demoBioAtomicLevel"));
    await session.step(21, "Then the table should have 6 rows", () => rowCount(page, 6));
    await session.step(22, "When user double-clicks on the \"cell 1 of molfile(HELM)\" area of grid", () => doubleClickArea(page, "cell 1 of molfile(HELM)", el("grid")));
    await session.step(23, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(24, "And the \"atoms\" reading of crux sketcher widget should be 339", () => readingIs(page, "atoms", el("crux sketcher widget"), 339));
    await session.step(25, "And crux sketcher widget should have an \"atom 338\" area", () => hasArea(page, el("crux sketcher widget"), "atom 338"));
    await session.step(26, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(27, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(28, "Then the \"atoms\" reading of crux sketcher widget should be 340", () => readingIs(page, "atoms", el("crux sketcher widget"), 340));
    await session.step(29, "When user clicks on crux undo button", () => clickOn(page, el("crux undo button")));
    await session.step(30, "Then the \"atoms\" reading of crux sketcher widget should be 339", () => readingIs(page, "atoms", el("crux sketcher widget"), 339));
    await session.step(31, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(32, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
  });
  test("R1 and R2 drawn in Crux in the monomer manager fill its R-groups grid with R1 and R2", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(36, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(37, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(38, "When user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
    await session.step(39, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
    await session.step(43, "When user clicks on \"Options\" icon in monomer sketcher", () => clickOn(page, el("\"Options\" icon in monomer sketcher")));
    await session.step(44, "And user picks \"OpenChemLib\" from the open menu", () => pickFromOpenMenu(page, "OpenChemLib"));
    await session.step(45, "And user clicks on \"Options\" icon in monomer sketcher", () => clickOn(page, el("\"Options\" icon in monomer sketcher")));
    await session.step(46, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(47, "And user clicks on \"Add New Monomer\" icon", () => clickOn(page, el("\"Add New Monomer\" icon")));
    await session.step(48, "Then the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(49, "And the \"atoms\" reading of crux sketcher widget should be 0", () => readingIs(page, "atoms", el("crux sketcher widget"), 0));
    await session.step(50, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(51, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(52, "And user clicks on the \"atom 1\" area of crux sketcher widget", () => clickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(53, "And user clicks on crux R-group tool", () => clickOn(page, el("crux R-group tool")));
    await session.step(54, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(55, "And user clicks on crux R1 button", () => clickOn(page, el("crux R1 button")));
    await session.step(56, "And user clicks on crux R-Group OK button", () => clickOn(page, el("crux R-Group OK button")));
    await session.step(57, "And user clicks on the \"atom 2\" area of crux sketcher widget", () => clickArea(page, "atom 2", el("crux sketcher widget")));
    await session.step(58, "And user clicks on crux R2 button", () => clickOn(page, el("crux R2 button")));
    await session.step(59, "And user clicks on crux R-Group OK button", () => clickOn(page, el("crux R-Group OK button")));
    await session.step(60, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]C[*:2]\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]C[*:2]"));
    await session.step(61, "When user clicks on \"R-groups\" tab", () => clickOn(page, el("\"R-groups\" tab")));
    await session.step(62, "Then the R-groups grid of the monomer form should list \"R1, R2\"", () => rgroupLabels(page, "R1, R2"));
  });
});
