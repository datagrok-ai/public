/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/fix-compound-form.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {compoundCount, compoundsNamed, noCompound, openEmsView} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: fix what the compound form names", () => {
  const session = feature(test, "features/guides/fix-compound-form.feature", import.meta.url);
  test("Fix a mistyped CAS and a missing structure, then register the compound", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page), undefined, "Sign in");
    await session.step(12, "And simple mode is off", () => simpleModeOff(page), undefined, "Simple mode is off");
    await session.step(14, "And no compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" is registered, now and at feature end", () => noCompound(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "Make sure the compound is not registered yet");
    await session.step(16, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"), undefined, "Open the Compounds view of MolTrack");
    enter(page, "MolTrack EMS");
    await session.step(18, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Click New to start a compound");
    await session.step(20, "And user types \"Triacontanoic acid\" into \"Common Name\" text input", () => typeInto(page, "Triacontanoic acid", el("\"Common Name\" text input")), undefined, "Type its common name");
    await session.step(22, "And user types \"506-50\" into \"CAS\" text input", () => typeInto(page, "506-50", el("\"CAS\" text input")), undefined, "Now a mistyped CAS number");
    await session.step(24, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save");
    await session.step(26, "Then save problem should contain text \"CAS must look like 50-78-2\"", () => shouldContainText(page, el("save problem"), "CAS must look like 50-78-2"), undefined, "Nothing is saved, and the status bar says why");
    await session.step(28, "And error of \"CAS\" text input should be visible", () => shouldBe(page, el("error of \"CAS\" text input"), "visible"), undefined, "The CAS field shows what is wrong");
    await session.step(30, "And error of \"Structure\" dart input should be visible", () => shouldBe(page, el("error of \"Structure\" dart input"), "visible"), undefined, "The empty Structure is flagged too, though it was never touched");
    await session.step(32, "And 0 compounds named \"Triacontanoic acid\" should be registered on the server", () => compoundsNamed(page, 0, "Triacontanoic acid"), undefined, "Nothing reached the server");
    await session.step(34, "When user types \"506-50-3\" into \"CAS\" text input", () => typeInto(page, "506-50-3", el("\"CAS\" text input")), undefined, "Fix the CAS number");
    await session.step(36, "Then error of \"CAS\" text input should be hidden", () => shouldBe(page, el("error of \"CAS\" text input"), "hidden"), undefined, "Its message goes away");
    await session.step(38, "When user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")), undefined, "Click the structure box to open the sketcher");
    await session.step(40, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O", el("molecule input of sketcher dialog")), undefined, "Paste the SMILES");
    await session.step(42, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
    await session.step(44, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept the structure");
    await session.step(46, "Then compound identity should contain text \"C30H60O2\"", () => shouldContainText(page, el("compound identity"), "C30H60O2"), undefined, "The formula is worked out right away");
    await session.step(48, "And error of \"Structure\" dart input should be hidden", () => shouldBe(page, el("error of \"Structure\" dart input"), "hidden"), undefined, "The Structure is no longer flagged");
    await session.step(50, "When user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save again");
    await session.step(52, "Then notification should contain text \"Compound saved\"", () => shouldContainText(page, el("notification"), "Compound saved"), undefined, "The compound is saved");
    await session.step(54, "And breadcrumbs should not contain text \"New compound\"", () => shouldNotContainText(page, el("breadcrumbs"), "New compound"), undefined, "It now has its Grok ID");
    await session.step(56, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "And it is on the server");
  });
});
