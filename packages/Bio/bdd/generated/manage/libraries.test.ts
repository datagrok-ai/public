/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/manage/libraries.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.menu.manage.monomer-libraries, bio.op.load_monomer_library, bio.op.save_monomer_library]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {libraryFileConforms} from '../../bindings/library-files.js';
import {allLibrariesSelected, chooseLibraryStorage, knownMonomer, loadedFrom, noLibraryOnServer, noSuchLibrary, notLoadedFrom, unknownMonomer} from '../../bindings/monomer-libs.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, pressKey, shouldBe, shouldContainText, uncheck, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {customFired, listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {call} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {closeCurrentView, openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Managing monomer libraries", () => {
  const session = feature(test, "features/manage/libraries.feature", import.meta.url);
  test("Managing monomer libraries", {tag: ["@journey", "@serial", "@realizes:bio.menu.manage.monomer-libraries", "@realizes:bio.op.load_monomer_library", "@realizes:bio.op.save_monomer_library"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(23, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(24, "And all monomer libraries are selected", () => allLibrariesSelected(page));
    await session.step(25, "And no \"bdd-test-lib.json\" monomer library is on the server", () => noLibraryOnServer(page, "bdd-test-lib.json"));
    await session.step(26, "When user picks \"Bio > Manage > Monomer Libraries\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomer Libraries"));
    await session.step(27, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(28, "And the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
    await session.step(29, "And \"Manage Duplicate Monomer Symbols\" heading should be visible", () => shouldBe(page, el("\"Manage Duplicate Monomer Symbols\" heading"), "visible"));
    await session.step(30, "And \"HELMCoreLibrary.json\" checkbox should be checked", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "checked"));
    await session.step(31, "And Search input should be visible", () => shouldBe(page, el("Search input"), "visible"));
    await session.step(32, "And Add button should be visible", () => shouldBe(page, el("Add button"), "visible"));
    await session.step(33, "And Merge button should be visible", () => shouldBe(page, el("Merge button"), "visible"));
    await run.scenario("The shipped library file is a HELM monomer library", async () => {
      await session.step(36, "Then the \"HELMCoreLibrary.json\" monomer library file should list monomers with a symbol and a structure each", () => libraryFileConforms(page, "HELMCoreLibrary.json"));
      await session.step(37, "And the monomer library should be loaded from \"HELMCoreLibrary.json\"", () => loadedFrom(page, "HELMCoreLibrary.json"));
    });
    await run.scenario("Unchecking a library reloads the monomer library without it, checking it back", async () => {
      await session.step(40, "Given user listens for \"bio-monomer-lib-loaded\" custom event", () => listenCustom(page, "bio-monomer-lib-loaded"));
      await session.step(41, "When user unchecks \"HELMCoreLibrary.json\" checkbox", () => uncheck(page, el("\"HELMCoreLibrary.json\" checkbox")));
      await session.step(42, "Then the \"bio-monomer-lib-loaded\" custom event should have fired", () => customFired(page, "bio-monomer-lib-loaded"));
      await session.step(43, "And the monomer library should not be loaded from \"HELMCoreLibrary.json\"", () => notLoadedFrom(page, "HELMCoreLibrary.json"));
      await session.step(44, "And \"HELMCoreLibrary.json\" checkbox should be unchecked", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "unchecked"));
      await session.step(45, "And \"polytool-lib.json\" checkbox should be checked", () => shouldBe(page, el("\"polytool-lib.json\" checkbox"), "checked"));
      await session.step(46, "And the monomer library should be loaded from \"polytool-lib.json\"", () => loadedFrom(page, "polytool-lib.json"));
      await session.step(47, "When user checks \"HELMCoreLibrary.json\" checkbox", () => check(page, el("\"HELMCoreLibrary.json\" checkbox")));
      await session.step(48, "Then the \"bio-monomer-lib-loaded\" custom event should have fired", () => customFired(page, "bio-monomer-lib-loaded"));
      await session.step(49, "And the monomer library should be loaded from \"HELMCoreLibrary.json\"", () => loadedFrom(page, "HELMCoreLibrary.json"));
      await session.step(50, "And \"A\" should be a known \"PEPTIDE\" monomer", () => knownMonomer(page, "A", "PEPTIDE"));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Add uploads a library file and its monomers become known", async () => {
      await session.step(54, "Given user listens for \"bio-monomer-lib-loaded\" custom event", () => listenCustom(page, "bio-monomer-lib-loaded"));
      await session.step(55, "And \"BDD\" should not be a known \"PEPTIDE\" monomer", () => unknownMonomer(page, "BDD", "PEPTIDE"));
      await session.step(56, "When user uploads \"fixtures/bdd-test-lib.json\" through Add button", () => uploadThrough(page, "fixtures/bdd-test-lib.json", el("Add button")));
      await session.step(57, "And user chooses \"Files\" storage for the uploaded monomer library", () => chooseLibraryStorage(page, "Files"));
      await session.step(58, "Then \"bdd-test-lib.json\" checkbox should become visible", () => shouldBe(page, el("\"bdd-test-lib.json\" checkbox"), "visible"));
      await session.step(59, "And \"bdd-test-lib.json\" checkbox should be checked", () => shouldBe(page, el("\"bdd-test-lib.json\" checkbox"), "checked"));
      await session.step(60, "And the \"bio-monomer-lib-loaded\" custom event should have fired", () => customFired(page, "bio-monomer-lib-loaded"));
      await session.step(61, "And the monomer library should be loaded from \"bdd-test-lib.json\"", () => loadedFrom(page, "bdd-test-lib.json"));
      await session.step(62, "And \"BDD\" should be a known \"PEPTIDE\" monomer", () => knownMonomer(page, "BDD", "PEPTIDE"));
      await session.step(63, "And an info balloon containing \"Added bdd-test-lib.json HELM library\" should have been shown", () => infoBalloonText(page, "Added bdd-test-lib.json HELM library"));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The uploaded library is still listed when the manager opens again", async () => {
      await session.step(67, "When user closes the current view", () => closeCurrentView(page));
      await session.step(68, "And user picks \"Bio > Manage > Monomer Libraries\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomer Libraries"));
      await session.step(69, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
      await session.step(70, "And \"bdd-test-lib.json\" checkbox should be visible", () => shouldBe(page, el("\"bdd-test-lib.json\" checkbox"), "visible"));
      await session.step(71, "And \"bdd-test-lib.json\" checkbox should be checked", () => shouldBe(page, el("\"bdd-test-lib.json\" checkbox"), "checked"));
      await session.step(72, "And \"HELMCoreLibrary.json\" checkbox should be checked", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "checked"));
    });
    await run.scenario("Delete unloads the library and removes its file after confirmation", async () => {
      await session.step(75, "Given user listens for \"bio-monomer-lib-loaded\" custom event", () => listenCustom(page, "bio-monomer-lib-loaded"));
      await session.step(76, "When user clicks on \"Delete\" icon in \"bdd-test-lib.json\" checkbox", () => clickOn(page, el("\"Delete\" icon in \"bdd-test-lib.json\" checkbox")));
      await session.step(77, "Then \"Warning\" dialog should be visible", () => shouldBe(page, el("\"Warning\" dialog"), "visible"));
      await session.step(78, "And \"Warning\" dialog should contain text \"bdd-test-lib.json\"", () => shouldContainText(page, el("\"Warning\" dialog"), "bdd-test-lib.json"));
      await session.step(79, "When user clicks on OK button in \"Warning\" dialog", () => clickOn(page, el("OK button in \"Warning\" dialog")));
      await session.step(80, "Then \"bdd-test-lib.json\" checkbox should become absent", () => shouldBe(page, el("\"bdd-test-lib.json\" checkbox"), "absent"));
      await session.step(81, "And there should be no \"bdd-test-lib.json\" monomer library on the server", () => noSuchLibrary(page, "bdd-test-lib.json"));
      await session.step(82, "And the \"bio-monomer-lib-loaded\" custom event should have fired", () => customFired(page, "bio-monomer-lib-loaded"));
      await session.step(83, "And the monomer library should not be loaded from \"bdd-test-lib.json\"", () => notLoadedFrom(page, "bdd-test-lib.json"));
      await session.step(84, "And \"BDD\" should not be a known \"PEPTIDE\" monomer", () => unknownMonomer(page, "BDD", "PEPTIDE"));
      await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The dialog entry lists the same libraries as the view", async () => {
      await session.step(88, "When user closes the current view", () => closeCurrentView(page));
      await session.step(89, "And user calls \"Bio:manageMonomerLibraries\" function", () => call(page, "Bio:manageMonomerLibraries"));
      await session.step(90, "Then \"Manage monomer libraries\" dialog should be visible", () => shouldBe(page, el("\"Manage monomer libraries\" dialog"), "visible"));
      await session.step(91, "And \"HELMCoreLibrary.json\" checkbox in \"Manage monomer libraries\" dialog should be checked", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox in \"Manage monomer libraries\" dialog"), "checked"));
      await session.step(92, "And \"polytool-lib.json\" checkbox in \"Manage monomer libraries\" dialog should be checked", () => shouldBe(page, el("\"polytool-lib.json\" checkbox in \"Manage monomer libraries\" dialog"), "checked"));
      await session.step(93, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(94, "Then \"Manage monomer libraries\" dialog should be absent", () => shouldBe(page, el("\"Manage monomer libraries\" dialog"), "absent"));
    });
    await run.scenario("After the dialog, the view lists the libraries again", async () => {
      await session.step(97, "When user picks \"Bio > Manage > Monomer Libraries\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomer Libraries"));
      await session.step(98, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
      await session.step(99, "And \"HELMCoreLibrary.json\" checkbox should be visible", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "visible"));
      await session.step(100, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
