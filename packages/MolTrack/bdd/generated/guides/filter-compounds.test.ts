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

test.describe("MolTrack: filter compounds by structure and properties", () => {
  const session = feature(test, "features/guides/filter-compounds.feature", import.meta.url);
  test("Phenols or anilines, without chlorine", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the MolTrack demo data is on the server", () => demoData(page));
    await session.step(13, "And user opens the \"Compounds\" view of MolTrack", () => openEmsView(page, "Compounds"));
    enter(page, "MolTrack");
    await session.step(14, "When user remembers the text of compound count", () => rememberText(page, el("compound count")));
    await session.step(16, "And user clicks on advanced mode in filters pane", () => clickOn(page, el("advanced mode in filters pane")), undefined, "The filters are on the left; Advanced nests groups of conditions");
    await session.step(18, "And user clicks on add group in filters pane", () => clickOn(page, el("add group in filters pane")), undefined, "A group inside: either of two substructures");
    await session.step(19, "And user clicks on last add condition in filters pane", () => clickOn(page, el("last add condition in filters pane")));
    await session.step(20, "And user selects \"Parent Structure\" in first condition property", () => selectIn(page, "Parent Structure", el("first condition property")));
    await session.step(21, "And user selects \"has substructure\" in first condition operator", () => selectIn(page, "has substructure", el("first condition operator")));
    await session.step(23, "And user clicks on molecule canvas in first condition value", () => clickOn(page, el("molecule canvas in first condition value")), undefined, "The first: a phenol, drawn in the sketcher");
    await session.step(24, "And user types \"c1ccccc1O\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1O", el("molecule input of sketcher dialog")));
    await session.step(25, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(26, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(28, "And user clicks on last add condition in filters pane", () => clickOn(page, el("last add condition in filters pane")), undefined, "The second, OR the first: an aniline");
    await session.step(29, "And user selects \"Parent Structure\" in second condition property", () => selectIn(page, "Parent Structure", el("second condition property")));
    await session.step(30, "And user selects \"has substructure\" in second condition operator", () => selectIn(page, "has substructure", el("second condition operator")));
    await session.step(31, "And user clicks on molecule canvas in second condition value", () => clickOn(page, el("molecule canvas in second condition value")));
    await session.step(32, "And user types \"c1ccccc1N\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1N", el("molecule input of sketcher dialog")));
    await session.step(33, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(34, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(35, "Then compound count should not have the remembered text", () => lostRememberedText(page, el("compound count")));
    await session.step(36, "And compound list should contain text \"Triclosan\"", () => shouldContainText(page, el("compound list"), "Triclosan"));
    await session.step(38, "When user clicks on first add condition in filters pane", () => clickOn(page, el("first add condition in filters pane")), undefined, "AND, outside the group: no chlorine in the formula");
    await session.step(39, "And user selects \"Formula\" in third condition property", () => selectIn(page, "Formula", el("third condition property")));
    await session.step(40, "And user selects \"does not contain\" in third condition operator", () => selectIn(page, "does not contain", el("third condition operator")));
    await session.step(41, "And user types \"Cl\" into third condition value", () => typeInto(page, "Cl", el("third condition value")));
    await session.step(42, "And user presses Tab in third condition value", () => pressKeyIn(page, "Tab", el("third condition value")));
    await session.step(44, "Then compound list should contain text \"Bisphenol A\"", () => shouldContainText(page, el("compound list"), "Bisphenol A"), undefined, "Phenols and anilines without chlorine: triclosan is gone");
    await session.step(45, "And compound list should not contain text \"Triclosan\"", () => shouldNotContainText(page, el("compound list"), "Triclosan"));
    await session.step(47, "When user clicks on third condition remove", () => clickOn(page, el("third condition remove")), undefined, "Remove the chlorine condition: triclosan is back");
    await session.step(48, "Then compound list should contain text \"Triclosan\"", () => shouldContainText(page, el("compound list"), "Triclosan"));
  });
});
