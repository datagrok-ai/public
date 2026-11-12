/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/atom-properties.feature
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
import {clickOn, selectIn, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Charges, isotopes and radicals in Crux", () => {
  const session = feature(test, "features/guides/crux/atom-properties.feature", import.meta.url);
  test("Make propanoate with the charge tool, label its methyl 13C and make the middle carbon a radical", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And the Crux sketcher is open on \"CCC(=O)O\"", () => cruxOpenOn(page, "CCC(=O)O"));
    await session.step(15, "When user clicks on Crux charge minus tool", () => clickOn(page, el("Crux charge minus tool")), undefined, "Pick the charge minus tool");
    await session.step(17, "And user clicks on the \"atom 4\" area of Crux sketcher widget", () => clickArea(page, "atom 4", el("Crux sketcher widget")), undefined, "Click the OH oxygen: it becomes O⁻");
    await session.step(19, "Then the sketcher in sketcher dialog should hold the molecule \"CCC(=O)[O-]\"", () => sketcherHolds(page, el("sketcher dialog"), "CCC(=O)[O-]"), undefined, "Propanoate");
    await session.step(21, "When user opens the Crux context menu on the \"atom 0\" area", () => cruxMenu(page, "atom 0"), undefined, "Right-click the methyl carbon");
    await session.step(23, "And user clicks on Crux Atom Properties item", () => clickOn(page, el("Crux Atom Properties item")), undefined, "Choose Atom Properties…");
    await session.step(25, "And user types \"13\" into Crux isotope field", () => typeInto(page, "13", el("Crux isotope field")), undefined, "Set the isotope to 13");
    await session.step(27, "And user clicks on Crux Atom Properties Apply button", () => clickOn(page, el("Crux Atom Properties Apply button")), undefined, "Apply");
    await session.step(29, "Then the sketcher in sketcher dialog should hold the molecule \"[13CH3]CC(=O)[O-]\"", () => sketcherHolds(page, el("sketcher dialog"), "[13CH3]CC(=O)[O-]"), undefined, "The methyl is now ¹³C");
    await session.step(31, "When user opens the Crux context menu on the \"atom 1\" area", () => cruxMenu(page, "atom 1"), undefined, "Right-click the middle carbon");
    await session.step(33, "And user clicks on Crux Atom Properties item", () => clickOn(page, el("Crux Atom Properties item")), undefined, "Choose Atom Properties…");
    await session.step(35, "And user selects \"Monoradical\" in Crux radical list", () => selectIn(page, "Monoradical", el("Crux radical list")), undefined, "Make it a monoradical: one unpaired electron");
    await session.step(37, "And user clicks on Crux Atom Properties Apply button", () => clickOn(page, el("Crux Atom Properties Apply button")), undefined, "Apply");
    await session.step(39, "Then the sketcher in sketcher dialog should hold the molecule \"[13CH3][CH]C(=O)[O-]\"", () => sketcherHolds(page, el("sketcher dialog"), "[13CH3][CH]C(=O)[O-]"), undefined, "The middle carbon now carries a radical");
  });
});
