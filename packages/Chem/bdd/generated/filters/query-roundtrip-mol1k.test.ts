/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/query-roundtrip-mol1k.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {filterPassesMatching, ketcherQueryFinds, readingFindsAs} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A query drawn in the substructure filter on mol1K: its rows, the sketcher reopened, and Ketcher", () => {
  const session = feature(test, "features/filters/query-roundtrip-mol1k.feature", import.meta.url);
  test("Toluene drawn in Crux, its methyl bond then marked aromatic, filters by the aromatic bond, and the sketcher reopened, Ketcher and Crux again hold that bond", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "Given user opens mol1K.sdf dataset", () => openDataset(page, ds("mol1K.sdf")));
    await session.step(21, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(22, "And user clicks on \"Sketch\" text in \"molecule\" filter card", () => clickOn(page, el("\"Sketch\" text in \"molecule\" filter card")));
    await session.step(23, "And user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(24, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(25, "And user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(26, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(27, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"Cc1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "Cc1ccccc1"));
    await session.step(28, "And the filter should pass exactly the molecules of \"molecule\" column containing \"[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => filterPassesMatching(page, "molecule", "[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"));
    await session.step(29, "When user clicks on crux aromatic bond tool", () => clickOn(page, el("crux aromatic bond tool")));
    await session.step(30, "And user clicks on the \"bond 6\" area of crux sketcher widget", () => clickArea(page, "bond 6", el("crux sketcher widget")));
    await session.step(31, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(32, "Then the filter should pass exactly the molecules of \"molecule\" column containing \"[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => filterPassesMatching(page, "molecule", "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"));
    await session.step(33, "And 293 rows should pass the filter", () => filterPasses(page, 293));
    await session.step(34, "When user clicks on sketcher thumbnail in \"molecule\" filter card", () => clickOn(page, el("sketcher thumbnail in \"molecule\" filter card")));
    await session.step(35, "Then the \"smarts\" reading of crux sketcher widget should find exactly the molecules of \"molecule\" column containing \"[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => readingFindsAs(page, "smarts", el("crux sketcher widget"), "molecule", "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"));
    await session.step(36, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(37, "And user picks \"Ketcher\" from the open menu", () => pickFromOpenMenu(page, "Ketcher"));
    await session.step(38, "Then the query Ketcher shows in sketcher dialog should find exactly the molecules of \"molecule\" column containing \"[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => ketcherQueryFinds(page, "molecule", "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"));
    await session.step(39, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(40, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(41, "Then the \"smarts\" reading of crux sketcher widget should find exactly the molecules of \"molecule\" column containing \"[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => readingFindsAs(page, "smarts", el("crux sketcher widget"), "molecule", "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"));
  });
});
