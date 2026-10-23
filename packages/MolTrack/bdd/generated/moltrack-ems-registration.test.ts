/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/moltrack-ems-registration.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {compoundCount, compoundsNamed, hasRememberedText, namesRemembered, noCompound, openEms, rememberCompound, rememberText, serverRefusesDuplicate} from '../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: registering a compound", () => {
  const session = feature(test, "features/moltrack-ems-registration.feature", import.meta.url);
  test("MolTrack EMS: registering a compound", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page), undefined, "Sign in");
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"), undefined, "Chem is installed");
    await session.step(16, "And no compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" is registered, now and at feature end", () => noCompound(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "Make sure the test compound is not registered yet");
    await session.step(18, "And user opens the MolTrack EMS app", () => openEms(page), undefined, "Open MolTrack");
    enter(page, "MolTrack EMS");
    await run.scenario("Save names every field an incomplete compound still needs, and writes nothing", async () => {
      await session.step(22, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Click New to start a compound");
      await session.step(24, "Then context panel should not contain text \"DG-\"", () => shouldNotContainText(page, el("context panel"), "DG-"), undefined, "A new compound page opens");
      await session.step(26, "When user types \"Octacosanoic acid\" into \"Common Name\" text input", () => typeInto(page, "Octacosanoic acid", el("\"Common Name\" text input")), undefined, "Type a common name");
      await session.step(28, "And user types \"123\" into \"CAS\" text input", () => typeInto(page, "123", el("\"CAS\" text input")), undefined, "Type a CAS number with a mistake");
      await session.step(30, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save");
      await session.step(32, "Then save problem should contain text \"CAS must look like 50-78-2\"", () => shouldContainText(page, el("save problem"), "CAS must look like 50-78-2"), undefined, "Nothing is saved, and the status bar says why");
      await session.step(34, "And error of \"CAS\" text input should be visible", () => shouldBe(page, el("error of \"CAS\" text input"), "visible"), undefined, "The CAS field shows what is wrong");
      await session.step(36, "And error of \"Structure\" dart input should be visible", () => shouldBe(page, el("error of \"Structure\" dart input"), "visible"), undefined, "The empty Structure is flagged too");
      await session.step(38, "And error of \"Structure\" dart input should contain text \"Value can't be empty\"", () => shouldContainText(page, el("error of \"Structure\" dart input"), "Value can't be empty"), undefined, "Its message reads value can't be empty");
      await session.step(40, "When user types \"506-48-9\" into \"CAS\" text input", () => typeInto(page, "506-48-9", el("\"CAS\" text input")), undefined, "Fix the CAS number");
      await session.step(42, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save again");
      await session.step(44, "Then save problem should contain text \"Structure is required\"", () => shouldContainText(page, el("save problem"), "Structure is required"), undefined, "Now the status bar says the structure is required");
      await session.step(46, "And 0 compounds named \"Octacosanoic acid\" should be registered on the server", () => compoundsNamed(page, 0, "Octacosanoic acid"), undefined, "Nothing reached the server");
      await session.step(48, "When user clicks on \"Discard\" button in app ribbon", () => clickOn(page, el("\"Discard\" button in app ribbon")), undefined, "Discard the draft");
      await session.step(50, "Then domain list should be visible", () => shouldBe(page, el("domain list"), "visible"), undefined, "Back on the list");
    });
    await run.scenario("A new structure is registered through the form", async () => {
      await session.step(54, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Click New again");
      await session.step(56, "Then breadcrumbs should contain text \"New compound\"", () => shouldContainText(page, el("breadcrumbs"), "New compound"), undefined, "The page is a new compound");
      await session.step(58, "When user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")), undefined, "Open the sketcher");
      await session.step(60, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O", el("molecule input of sketcher dialog")), undefined, "Paste the SMILES");
      await session.step(62, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
      await session.step(64, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept the structure");
      await session.step(66, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"), undefined, "The sketcher is closed");
      await session.step(68, "And compound identity should contain text \"C28H56O2\"", () => shouldContainText(page, el("compound identity"), "C28H56O2"), undefined, "The formula is worked out right away");
      await session.step(70, "When user types \"Octacosanoic acid\" into \"Common Name\" text input", () => typeInto(page, "Octacosanoic acid", el("\"Common Name\" text input")), undefined, "Type its common name");
      await session.step(72, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save");
      await session.step(74, "Then notification should contain text \"Compound saved\"", () => shouldContainText(page, el("notification"), "Compound saved"), undefined, "The compound is saved");
      await session.step(76, "And breadcrumbs should not contain text \"New compound\"", () => shouldNotContainText(page, el("breadcrumbs"), "New compound"), undefined, "It now has its Grok ID");
      await session.step(78, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "And it is on the server");
      await session.step(80, "When user remembers the compound the breadcrumbs name", () => rememberCompound(page), undefined, "Remember its Grok ID");
    });
    await run.scenario("Its sodium salt is refused as the compound already registered", async () => {
      await session.step(84, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Now its sodium salt: click New");
      await session.step(86, "And user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")), undefined, "Open the sketcher");
      await session.step(88, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]", el("molecule input of sketcher dialog")), undefined, "Paste the salt");
      await session.step(90, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
      await session.step(92, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept the structure");
      await session.step(94, "Then compound identity should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("compound identity"), "Already registered as DG-"), undefined, "The identity panel says it is already registered");
      await session.step(96, "And compound identity should name the remembered compound", () => namesRemembered(page, el("compound identity")), undefined, "It names the compound registered a moment ago");
      await session.step(98, "And error of \"Structure\" dart input should be visible", () => shouldBe(page, el("error of \"Structure\" dart input"), "visible"), undefined, "The Structure field shows the message too");
      await session.step(100, "And error of \"Structure\" dart input should name the remembered compound", () => namesRemembered(page, el("error of \"Structure\" dart input")), undefined, "And names that compound");
      await session.step(102, "When user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save");
      await session.step(104, "Then save problem should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("save problem"), "Already registered as DG-"), undefined, "Save is refused, and says why");
      await session.step(106, "When user clicks on \"Discard\" button in app ribbon", () => clickOn(page, el("\"Discard\" button in app ribbon")), undefined, "Discard the draft");
      await session.step(108, "Then domain list should be visible", () => shouldBe(page, el("domain list"), "visible"), undefined, "Back on the list");
      await session.step(110, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "Still only one such compound on the server");
      await session.step(112, "And the server should refuse another compound with the remembered compound's registration hash", () => serverRefusesDuplicate(page), undefined, "The server itself refuses a second compound with the same identity");
    });
    await run.scenario("A substructure condition narrows the list to the compound", async () => {
      await session.step(116, "When user remembers the text of compound count", () => rememberText(page, el("compound count")), undefined, "Note how many compounds the list shows");
      await session.step(118, "And user clicks on add condition in filters pane", () => clickOn(page, el("add condition in filters pane")), undefined, "Add a filter condition");
      await session.step(120, "And user selects \"Parent Structure\" in condition property", () => selectIn(page, "Parent Structure", el("condition property")), undefined, "Pick the parent structure");
      await session.step(122, "And user selects \"has substructure\" in condition operator", () => selectIn(page, "has substructure", el("condition operator")), undefined, "Choose has substructure");
      await session.step(124, "And user clicks on molecule canvas in condition value", () => clickOn(page, el("molecule canvas in condition value")), undefined, "Open the sketcher for the condition");
      await session.step(126, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", el("molecule input of sketcher dialog")), undefined, "Paste the substructure");
      await session.step(128, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
      await session.step(130, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept it");
      await session.step(132, "Then compound list should have 1 item in total", () => shouldHaveItemsInTotal(page, el("compound list"), 1), undefined, "One compound is left");
      await session.step(134, "And compound count should have text \"1 compound\"", () => shouldHaveText(page, el("compound count"), "1 compound"), undefined, "The count says one compound");
      await session.step(136, "And compound list should contain text \"Octacosanoic acid\"", () => shouldContainText(page, el("compound list"), "Octacosanoic acid"), undefined, "It is the acid");
      await session.step(138, "And compound list should name the remembered compound", () => namesRemembered(page, el("compound list")), undefined, "With the Grok ID it was given");
      await session.step(140, "When user clicks on condition remove", () => clickOn(page, el("condition remove")), undefined, "Remove the condition");
      await session.step(142, "Then compound count should have the remembered text", () => hasRememberedText(page, el("compound count")), undefined, "All compounds are back");
    });
    run.finish();
  });
});
