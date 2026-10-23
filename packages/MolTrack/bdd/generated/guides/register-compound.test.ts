/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/register-compound.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {compoundCount, namesRemembered, noCompound, openEmsView, rememberCompound} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register a compound", () => {
  const session = feature(test, "features/guides/register-compound.feature", import.meta.url);
  test("Register an acid, then see its sodium salt refused as the same compound", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page), undefined, "Sign in");
    await session.step(14, "And simple mode is off", () => simpleModeOff(page), undefined, "Simple mode is off");
    await session.step(16, "And no compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" is registered, now and at feature end", () => noCompound(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "Make sure the compound is not registered yet");
    await session.step(18, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"), undefined, "Open the Compounds view of MolTrack");
    enter(page, "MolTrack EMS");
    await session.step(20, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Click New to start a compound");
    await session.step(22, "And user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")), undefined, "Click the structure box to open the sketcher");
    await session.step(24, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O", el("molecule input of sketcher dialog")), undefined, "Paste the SMILES");
    await session.step(26, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
    await session.step(28, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept the structure");
    await session.step(30, "Then compound identity should contain text \"C29H58O2\"", () => shouldContainText(page, el("compound identity"), "C29H58O2"), undefined, "The formula and the identity are worked out before Save");
    await session.step(32, "When user types \"Nonacosanoic acid\" into \"Common Name\" text input", () => typeInto(page, "Nonacosanoic acid", el("\"Common Name\" text input")), undefined, "Give it a common name");
    await session.step(34, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save");
    await session.step(36, "Then notification should contain text \"Compound saved\"", () => shouldContainText(page, el("notification"), "Compound saved"), undefined, "The compound is saved");
    await session.step(38, "And breadcrumbs should not contain text \"New compound\"", () => shouldNotContainText(page, el("breadcrumbs"), "New compound"), undefined, "It now has its Grok ID");
    await session.step(40, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "And it is on the server");
    await session.step(42, "When user remembers the compound the breadcrumbs name", () => rememberCompound(page), undefined, "Remember its Grok ID");
    await session.step(44, "And user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Now try its sodium salt: click New");
    await session.step(46, "And user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")), undefined, "Open the sketcher");
    await session.step(48, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]", el("molecule input of sketcher dialog")), undefined, "Paste the salt");
    await session.step(50, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
    await session.step(52, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept the structure");
    await session.step(54, "Then error of \"Structure\" dart input should name the remembered compound", () => namesRemembered(page, el("error of \"Structure\" dart input")), undefined, "The Structure field says it is the compound registered a moment ago");
    await session.step(56, "And compound identity should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("compound identity"), "Already registered as DG-"), undefined, "The identity panel says the same");
    await session.step(58, "When user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")), undefined, "Press Save anyway");
    await session.step(60, "Then save problem should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("save problem"), "Already registered as DG-"), undefined, "Save is refused, and the status bar says why");
    await session.step(62, "When user clicks on \"Discard\" button in app ribbon", () => clickOn(page, el("\"Discard\" button in app ribbon")), undefined, "Discard the draft");
    await session.step(64, "Then 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"), undefined, "Still only one such compound on the server");
  });
});
