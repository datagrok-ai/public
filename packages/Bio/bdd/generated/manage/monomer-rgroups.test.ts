/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/manage/monomer-rgroups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.menu.manage.monomers]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {monomerSketcherReady} from '../../bindings/library-files.js';
import {allLibrariesSelected, libraryOnServer} from '../../bindings/monomer-libs.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveError, shouldHaveNoError, shouldListItem} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, packageInstalled, sketcherIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, infoBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The R groups of a monomer in the Manage Monomers editor", () => {
  const session = feature(test, "features/manage/monomer-rgroups.feature", import.meta.url);
  test("The R groups of a monomer in the Manage Monomers editor", {tag: ["@journey", "@serial", "@realizes:bio.menu.manage.monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(22, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(23, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(24, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(25, "And all monomer libraries are selected", () => allLibrariesSelected(page));
    await session.step(26, "And the \"fixtures/bdd-rgroups-lib.json\" file of the project is the \"bdd-rgroups.json\" monomer library on the server", () => libraryOnServer(page, "fixtures/bdd-rgroups-lib.json", "bdd-rgroups.json"));
    await session.step(27, "When user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
    await session.step(28, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
    await session.step(29, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerSketcherReady(page));
    await session.step(30, "When user selects \"bdd-rgroups.json\" in \"Monomer Library\" input in monomer form", () => selectIn(page, "bdd-rgroups.json", el("\"Monomer Library\" input in monomer form")));
    await session.step(31, "Then the table should have 1 row", () => rowCount(page, 1));
    await run.scenario("Each R group shows the cap its definition means, however it is written", async () => {
      await session.step(34, "When user clicks on the \"cell 1 of Symbol\" area of grid", () => clickArea(page, "cell 1 of Symbol", el("grid")));
      await session.step(35, "And user clicks on \"Edit Monomer\" icon", () => clickOn(page, el("\"Edit Monomer\" icon")));
      await session.step(36, "And user clicks on \"R-groups\" tab", () => clickOn(page, el("\"R-groups\" tab")));
      await session.step(37, "Then \"R1\" input in monomer form should have value \"H\"", () => shouldHaveValue(page, el("\"R1\" input in monomer form"), "H"));
      await session.step(38, "And \"R2\" input in monomer form should have value \"OH\"", () => shouldHaveValue(page, el("\"R2\" input in monomer form"), "OH"));
      await session.step(39, "And \"R3\" input in monomer form should have value \"Azide\"", () => shouldHaveValue(page, el("\"R3\" input in monomer form"), "Azide"));
      await session.step(40, "And \"R4\" input in monomer form should list the item \"X\"", () => shouldListItem(page, el("\"R4\" input in monomer form"), "X"));
      await session.step(41, "And \"R4\" input in monomer form should have the error \"X is a placeholder, not a cap group — pick the cap R4 gets (usually H)\"", () => shouldHaveError(page, el("\"R4\" input in monomer form"), "X is a placeholder, not a cap group — pick the cap R4 gets (usually H)"));
      await session.step(42, "And \"R5\" input in monomer form should have the error \"Pick its cap group, or define it with the pencil\"", () => shouldHaveError(page, el("\"R5\" input in monomer form"), "Pick its cap group, or define it with the pencil"));
      await session.step(43, "And \"R2\" input in monomer form should contain the text \"hydroxyl · O[*:2]\"", () => shouldContainText(page, el("\"R2\" input in monomer form"), "hydroxyl · O[*:2]"));
      await session.step(44, "And \"Save Monomer\" button in monomer form should be disabled", () => shouldBe(page, el("\"Save Monomer\" button in monomer form"), "disabled"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pick and the pencil define the R groups, and Save writes them into the library", async () => {
      await session.step(48, "When user selects \"Br\" in \"R4\" input in monomer form", () => selectIn(page, "Br", el("\"R4\" input in monomer form")));
      await session.step(49, "Then \"R4\" input in monomer form should have no error", () => shouldHaveNoError(page, el("\"R4\" input in monomer form")));
      await session.step(50, "And \"R4\" input in monomer form should contain the text \"bromo · Br[*:4]\"", () => shouldContainText(page, el("\"R4\" input in monomer form"), "bromo · Br[*:4]"));
      await session.step(51, "When user clicks on \"Edit the R5 definition\" icon button", () => clickOn(page, el("\"Edit the R5 definition\" icon button")));
      await session.step(52, "And user enters \"[*:5]C(F)(F)F\" into \"Cap SMILES\" input in monomer form", () => enterInto(page, "[*:5]C(F)(F)F", el("\"Cap SMILES\" input in monomer form")));
      await session.step(53, "And user enters \"CF3\" into \"Cap name\" input in monomer form", () => enterInto(page, "CF3", el("\"Cap name\" input in monomer form")));
      await session.step(54, "And user enters \"R5-CF3\" into \"Alternate ID\" input in monomer form", () => enterInto(page, "R5-CF3", el("\"Alternate ID\" input in monomer form")));
      await session.step(55, "Then \"R5\" input in monomer form should have no error", () => shouldHaveNoError(page, el("\"R5\" input in monomer form")));
      await session.step(56, "And \"R5\" input in monomer form should list the item \"CF3\"", () => shouldListItem(page, el("\"R5\" input in monomer form"), "CF3"));
      await session.step(57, "And \"Save Monomer\" button in monomer form should be enabled", () => shouldBe(page, el("\"Save Monomer\" button in monomer form"), "enabled"));
      await session.step(58, "When user clicks on \"Save Monomer\" button in monomer form", () => clickOn(page, el("\"Save Monomer\" button in monomer form")));
      await session.step(59, "Then \"Save Monomer\" dialog should be visible", () => shouldBe(page, el("\"Save Monomer\" dialog"), "visible"));
      await session.step(60, "When user clicks on Save button in \"Save Monomer\" dialog", () => clickOn(page, el("Save button in \"Save Monomer\" dialog")));
      await session.step(61, "Then an info balloon containing \"Monomer bddRG was successfully saved in library bdd-rgroups.json\" should have been shown", () => infoBalloonText(page, "Monomer bddRG was successfully saved in library bdd-rgroups.json"));
      await session.step(62, "And the value of \"R1\" column in row 1 should be \"[*:1][H]\"", () => valueInRow(page, "R1", 1, "[*:1][H]"));
      await session.step(63, "And the value of \"R2\" column in row 1 should be \"O[*:2]\"", () => valueInRow(page, "R2", 1, "O[*:2]"));
      await session.step(64, "And the value of \"R3\" column in row 1 should be \"[N-]=[N+]=N[*:3]\"", () => valueInRow(page, "R3", 1, "[N-]=[N+]=N[*:3]"));
      await session.step(65, "And the value of \"R4\" column in row 1 should be \"Br[*:4]\"", () => valueInRow(page, "R4", 1, "Br[*:4]"));
      await session.step(66, "And the value of \"R5\" column in row 1 should be \"[*:5]C(F)(F)F\"", () => valueInRow(page, "R5", 1, "[*:5]C(F)(F)F"));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
