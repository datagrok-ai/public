/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/filter-compounds.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {demoData, lostRememberedText, openEmsView, rememberText} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, selectIn, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: filter compounds by structure and properties", () => {
  const session = feature(test, "features/guides/filter-compounds.feature", import.meta.url);
  test("Phenols or anilines, without chlorine", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page), undefined, "Sign in");
    await session.step(13, "And simple mode is off", () => simpleModeOff(page), undefined, "Simple mode is off");
    await session.step(15, "And the MolTrack EMS demo data is on the server", () => demoData(page), undefined, "The demo data is loaded");
    await session.step(17, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"), undefined, "Open the Compounds view of MolTrack");
    enter(page, "MolTrack EMS");
    await session.step(19, "When user remembers the text of compound count", () => rememberText(page, el("compound count")), undefined, "Note how many compounds there are");
    await session.step(21, "And user clicks on advanced mode in filters pane", () => clickOn(page, el("advanced mode in filters pane")), undefined, "Switch the filters to Advanced, to nest conditions");
    await session.step(23, "And user clicks on add group in filters pane", () => clickOn(page, el("add group in filters pane")), undefined, "Add a group: either of two substructures");
    await session.step(25, "And user clicks on last add condition in filters pane", () => clickOn(page, el("last add condition in filters pane")), undefined, "Add its first condition");
    await session.step(27, "And user selects \"Parent Structure\" in first condition property", () => selectIn(page, "Parent Structure", el("first condition property")), undefined, "On the parent structure");
    await session.step(29, "And user selects \"has substructure\" in first condition operator", () => selectIn(page, "has substructure", el("first condition operator")), undefined, "Has substructure");
    await session.step(31, "And user clicks on molecule canvas in first condition value", () => clickOn(page, el("molecule canvas in first condition value")), undefined, "Open the sketcher");
    await session.step(33, "And user types \"c1ccccc1O\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1O", el("molecule input of sketcher dialog")), undefined, "Type a phenol");
    await session.step(35, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
    await session.step(37, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept it");
    await session.step(39, "And user clicks on last add condition in filters pane", () => clickOn(page, el("last add condition in filters pane")), undefined, "Add the second condition to the group");
    await session.step(41, "And user selects \"Parent Structure\" in second condition property", () => selectIn(page, "Parent Structure", el("second condition property")), undefined, "Again on the parent structure");
    await session.step(43, "And user selects \"has substructure\" in second condition operator", () => selectIn(page, "has substructure", el("second condition operator")), undefined, "Has substructure");
    await session.step(45, "And user clicks on molecule canvas in second condition value", () => clickOn(page, el("molecule canvas in second condition value")), undefined, "Open the sketcher");
    await session.step(47, "And user types \"c1ccccc1N\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1N", el("molecule input of sketcher dialog")), undefined, "Type an aniline");
    await session.step(49, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")), undefined, "Press Enter to draw it");
    await session.step(51, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "Accept it: phenols or anilines");
    await session.step(53, "Then compound count should not have the remembered text", () => lostRememberedText(page, el("compound count")), undefined, "The list narrows");
    await session.step(55, "And compound list should contain text \"Triclosan\"", () => shouldContainText(page, el("compound list"), "Triclosan"), undefined, "Triclosan is among them");
    await session.step(57, "When user clicks on first add condition in filters pane", () => clickOn(page, el("first add condition in filters pane")), undefined, "Add a condition outside the group");
    await session.step(59, "And user selects \"Formula\" in third condition property", () => selectIn(page, "Formula", el("third condition property")), undefined, "On the formula");
    await session.step(61, "And user selects \"does not contain\" in third condition operator", () => selectIn(page, "does not contain", el("third condition operator")), undefined, "Does not contain");
    await session.step(63, "And user types \"Cl\" into third condition value", () => typeInto(page, "Cl", el("third condition value")), undefined, "Type Cl, for chlorine");
    await session.step(65, "And user presses Tab in third condition value", () => pressKeyIn(page, "Tab", el("third condition value")), undefined, "Press Tab to apply it");
    await session.step(67, "Then compound list should contain text \"Bisphenol A\"", () => shouldContainText(page, el("compound list"), "Bisphenol A"), undefined, "Bisphenol A stays");
    await session.step(69, "And compound list should not contain text \"Triclosan\"", () => shouldNotContainText(page, el("compound list"), "Triclosan"), undefined, "Triclosan, with chlorine, is gone");
    await session.step(71, "When user clicks on third condition remove", () => clickOn(page, el("third condition remove")), undefined, "Remove the chlorine condition");
    await session.step(73, "Then compound list should contain text \"Triclosan\"", () => shouldContainText(page, el("compound list"), "Triclosan"), undefined, "Triclosan is back");
  });
});
