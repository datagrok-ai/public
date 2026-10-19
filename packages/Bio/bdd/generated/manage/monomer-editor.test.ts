/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/manage/monomer-editor.feature
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
import {allLibrariesSelected, knownMonomer, libraryOnServer, unknownMonomer} from '../../bindings/monomer-libs.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, pressKeyIn, selectIn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, packageInstalled, sketcherIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, infoBalloonText, noBalloons, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Making and removing a monomer in the Manage Monomers view", () => {
  const session = feature(test, "features/manage/monomer-editor.feature", import.meta.url);
  test("Making and removing a monomer in the Manage Monomers view", {tag: ["@journey", "@serial", "@realizes:bio.menu.manage.monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(23, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(24, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(25, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(26, "And all monomer libraries are selected", () => allLibrariesSelected(page));
    await session.step(27, "And the \"fixtures/bdd-test-lib.json\" file of the project is the \"bdd-monomer-editor.json\" monomer library on the server", () => libraryOnServer(page, "fixtures/bdd-test-lib.json", "bdd-monomer-editor.json"));
    await session.step(28, "When user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
    await session.step(29, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
    await session.step(30, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerSketcherReady(page));
    await run.scenario("The Library input shows the monomers of the library picked", async () => {
      await session.step(33, "When user selects \"bdd-monomer-editor\" in \"Library\" input in monomer form", () => selectIn(page, "bdd-monomer-editor", el("\"Library\" input in monomer form")));
      await session.step(34, "Then the table should have 1 row", () => rowCount(page, 1));
      await session.step(35, "And the value of \"Symbol\" column in row 1 should be \"BDD\"", () => valueInRow(page, "Symbol", 1, "BDD"));
    });
    await run.scenario("A monomer drawn in the editor and saved joins the library", async () => {
      await session.step(38, "When user clicks on \"Add New Monomer\" icon", () => clickOn(page, el("\"Add New Monomer\" icon")));
      await session.step(39, "And user types \"CC(C)(C)[C@H](N[1*])C([2*])=O\" into molecule input of monomer sketcher", () => typeInto(page, "CC(C)(C)[C@H](N[1*])C([2*])=O", el("molecule input of monomer sketcher")));
      await session.step(40, "And user presses Enter in molecule input of monomer sketcher", () => pressKeyIn(page, "Enter", el("molecule input of monomer sketcher")));
      await session.step(41, "Then \"R1\" R-group chip in monomer form should contain the text \"R1 H\"", () => shouldContainText(page, el("\"R1\" R-group chip in monomer form"), "R1 H"));
      await session.step(42, "And \"R2\" R-group chip in monomer form should contain the text \"R2 OH\"", () => shouldContainText(page, el("\"R2\" R-group chip in monomer form"), "R2 OH"));
      await session.step(43, "When user enters \"BDDTLE\" into \"Symbol\" input in monomer form", () => enterInto(page, "BDDTLE", el("\"Symbol\" input in monomer form")));
      await session.step(44, "And user enters \"bdd tert-leucine\" into \"Name\" input in monomer form", () => enterInto(page, "bdd tert-leucine", el("\"Name\" input in monomer form")));
      await session.step(45, "And user clicks on \"Add Monomer\" button in monomer form", () => clickOn(page, el("\"Add Monomer\" button in monomer form")));
      await session.step(46, "Then an info balloon containing \"Monomer BDDTLE was successfully saved in library bdd-monomer-editor.json\" should have been shown", () => infoBalloonText(page, "Monomer BDDTLE was successfully saved in library bdd-monomer-editor.json"));
      await session.step(47, "And \"BDDTLE\" should be a known \"PEPTIDE\" monomer", () => knownMonomer(page, "BDDTLE", "PEPTIDE"));
      await session.step(48, "And the table should have 2 rows", () => rowCount(page, 2));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Delete removes the current monomer from the library after a confirmation", async () => {
      await session.step(53, "Given the value of \"Symbol\" column in row 2 should be \"BDDTLE\"", () => valueInRow(page, "Symbol", 2, "BDDTLE"));
      await session.step(54, "When user clicks on the \"cell 2 of Symbol\" area of grid", () => clickArea(page, "cell 2 of Symbol", el("grid")));
      await session.step(55, "Then the \"current row\" reading of grid should be 2", () => readingIs(page, "current row", el("grid"), 2));
      await session.step(56, "When user clicks on \"Delete\" icon", () => clickOn(page, el("\"Delete\" icon")));
      await session.step(57, "Then \"Remove Monomer\" dialog should be visible", () => shouldBe(page, el("\"Remove Monomer\" dialog"), "visible"));
      await session.step(58, "And \"Remove Monomer\" dialog should contain text \"BDDTLE\"", () => shouldContainText(page, el("\"Remove Monomer\" dialog"), "BDDTLE"));
      await session.step(59, "When user clicks on Remove button in \"Remove Monomer\" dialog", () => clickOn(page, el("Remove button in \"Remove Monomer\" dialog")));
      await session.step(60, "Then an info balloon containing \"Monomer BDDTLE was successfully removed from bdd-monomer-editor.json library\" should have been shown", () => infoBalloonText(page, "Monomer BDDTLE was successfully removed from bdd-monomer-editor.json library"));
      await session.step(61, "And \"BDDTLE\" should not be a known \"PEPTIDE\" monomer", () => unknownMonomer(page, "BDDTLE", "PEPTIDE"));
      await session.step(62, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
