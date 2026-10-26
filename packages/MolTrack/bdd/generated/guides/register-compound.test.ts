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
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And no compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" is registered, now and at feature end", () => noCompound(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"));
    await session.step(14, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"));
    enter(page, "MolTrack EMS");
    await session.step(16, "When user clicks on \"New compound\" button in app ribbon", () => clickOn(page, el("\"New compound\" button in app ribbon")), undefined, "Start a new compound");
    await session.step(18, "And user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")), undefined, "Draw the structure: click the empty picture to open the sketcher");
    await session.step(20, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O", el("molecule input of sketcher dialog")), undefined, "Type or paste a SMILES");
    await session.step(21, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(22, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(24, "Then compound identity should contain text \"C29H58O2\"", () => shouldContainText(page, el("compound identity"), "C29H58O2"), undefined, "The formula and the identity are computed before Save");
    await session.step(25, "When user types \"Nonacosanoic acid\" into \"Common Name\" text input", () => typeInto(page, "Nonacosanoic acid", el("\"Common Name\" text input")));
    await session.step(26, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(28, "Then notification should contain text \"Compound saved\"", () => shouldContainText(page, el("notification"), "Compound saved"), undefined, "Saved: the compound has its Grok ID");
    await session.step(29, "And breadcrumbs should not contain text \"New compound\"", () => shouldNotContainText(page, el("breadcrumbs"), "New compound"));
    await session.step(30, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"));
    await session.step(31, "When user remembers the compound the breadcrumbs name", () => rememberCompound(page));
    await session.step(33, "And user clicks on \"New compound\" button in app ribbon", () => clickOn(page, el("\"New compound\" button in app ribbon")), undefined, "Now its sodium salt");
    await session.step(34, "And user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")));
    await session.step(35, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]", el("molecule input of sketcher dialog")));
    await session.step(36, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(37, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(39, "Then compound identity should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("compound identity"), "Already registered as DG-"), undefined, "Refused: the salt is the compound registered a moment ago");
    await session.step(40, "And compound identity should name the remembered compound", () => namesRemembered(page, el("compound identity")));
    await session.step(41, "When user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(43, "Then save problem should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("save problem"), "Already registered as DG-"), undefined, "Save is refused, and says why");
    await session.step(44, "When user clicks on \"Discard\" button in app ribbon", () => clickOn(page, el("\"Discard\" button in app ribbon")));
    await session.step(45, "Then 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"));
  });
});
