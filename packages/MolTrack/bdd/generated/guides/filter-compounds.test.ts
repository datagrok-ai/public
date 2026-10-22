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
import {demoData, hasRememberedText, lostRememberedText, openEmsView, rememberText} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, pressKeyIn, selectIn, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: filter compounds by structure and properties", () => {
  const session = feature(test, "features/guides/filter-compounds.feature", import.meta.url);
  test("Phenols, then the compounds similar to bisphenol A", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(13, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"));
    enter(page, "MolTrack EMS");
    await session.step(14, "When user remembers the text of compound count", () => rememberText(page, el("compound count")));
    await session.step(16, "And user clicks on \"Show the filters\" button", () => clickOn(page, el("\"Show the filters\" button")), undefined, "Open the filters, and add a condition");
    await session.step(17, "And user clicks on add condition in filters pane", () => clickOn(page, el("add condition in filters pane")));
    await session.step(18, "And user selects \"Parent Structure\" in condition property", () => selectIn(page, "Parent Structure", el("condition property")));
    await session.step(19, "And user selects \"has substructure\" in condition operator", () => selectIn(page, "has substructure", el("condition operator")));
    await session.step(21, "And user clicks on molecule canvas in condition value", () => clickOn(page, el("molecule canvas in condition value")), undefined, "Draw the substructure: a phenol");
    await session.step(22, "And user types \"c1ccccc1O\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1O", el("molecule input of sketcher dialog")));
    await session.step(23, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(24, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(26, "Then compound count should not have the remembered text", () => lostRememberedText(page, el("compound count")), undefined, "Only the phenols are left");
    await session.step(27, "And compound list should contain text \"Triclosan\"", () => shouldContainText(page, el("compound list"), "Triclosan"));
    await session.step(29, "When user selects \"is similar to\" in condition operator", () => selectIn(page, "is similar to", el("condition operator")), undefined, "Now the compounds similar to bisphenol A, at 0.4 or more");
    await session.step(30, "And user clicks on molecule canvas in condition value", () => clickOn(page, el("molecule canvas in condition value")));
    await session.step(31, "And user types \"CC(C)(c1ccc(O)cc1)c1ccc(O)cc1\" into molecule input of sketcher dialog", () => typeInto(page, "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1", el("molecule input of sketcher dialog")));
    await session.step(32, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(33, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(34, "And user enters \"0.4\" into condition threshold", () => enterInto(page, "0.4", el("condition threshold")));
    await session.step(35, "Then compound list should contain text \"4-(1,1,3,3-Tetramethylbutyl)phenol\"", () => shouldContainText(page, el("compound list"), "4-(1,1,3,3-Tetramethylbutyl)phenol"));
    await session.step(37, "When user clicks on condition remove", () => clickOn(page, el("condition remove")), undefined, "Remove the condition: every compound is back");
    await session.step(38, "Then compound count should have the remembered text", () => hasRememberedText(page, el("compound count")));
  });
});
