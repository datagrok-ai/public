/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/abbreviations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/datasets.js';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cruxMenu, cruxOpenOn} from '../../../bindings/crux.js';
import {sketcherHolds} from '../../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Abbreviations in Crux", () => {
  const session = feature(test, "features/guides/crux/abbreviations.feature", import.meta.url);
  test("Type CO2Me on toluene's methyl, then expand the ester into its atoms", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"Cc1ccccc1\"", () => cruxOpenOn(page, "Cc1ccccc1"));
    await session.step(16, "When user double-clicks on the \"atom 0\" area of Crux sketcher widget", () => doubleClickArea(page, "atom 0", el("Crux sketcher widget")), undefined, "Double-click the methyl carbon to edit its label");
    await session.step(18, "And user types \"CO2Me\" into label editor of Crux sketcher widget", () => typeInto(page, "CO2Me", el("label editor of Crux sketcher widget")), undefined, "Type CO2Me");
    await session.step(20, "And user presses Enter in label editor of Crux sketcher widget", () => pressKeyIn(page, "Enter", el("label editor of Crux sketcher widget")), undefined, "Press Enter: a methyl ester, drawn as CO2Me");
    await session.step(22, "Then the sketcher in sketcher dialog should hold the molecule \"COC(=O)c1ccccc1\"", () => sketcherHolds(page, el("sketcher dialog"), "COC(=O)c1ccccc1"), undefined, "Methyl benzoate");
    await session.step(24, "When user opens the Crux context menu on the \"atom 0\" area", () => cruxMenu(page, "atom 0"), undefined, "Right-click the CO2Me label");
    await session.step(26, "And user clicks on Crux expand abbreviation item", () => clickOn(page, el("Crux expand abbreviation item")), undefined, "Choose Expand Abbreviation: the ester's atoms are drawn");
    await session.step(28, "Then the sketcher in sketcher dialog should hold the molecule \"COC(=O)c1ccccc1\"", () => sketcherHolds(page, el("sketcher dialog"), "COC(=O)c1ccccc1"), undefined, "The same molecule, drawn atom by atom");
  });
});
