/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-molecule-parameter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {scriptResult} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, doubleClickOn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, dialogCloses, packageInstalled, scriptOnServer, scriptsView, sketcherIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A script's molecule parameter, sketched in Crux", () => {
  const session = feature(test, "features/scripts/scripts-molecule-parameter.feature", import.meta.url);
  test("A pentane drawn in Crux for the Molecule parameter reaches the script as CCCCC", {tag: ["@serial", "@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And a script \"BddCruxMolecule{time}\" is on the server:", () => scriptOnServer(page, session.text("BddCruxMolecule{time}"), "//language: javascript\n//input: string mol {semType: Molecule}\n//output: string got\ngot = mol;"));
    await session.step(30, "Given user opens the Scripts view", () => scriptsView(page));
    await session.step(31, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(32, "And user types \"BddCruxMolecule{time}\" into gallery search", () => typeInto(page, session.text("BddCruxMolecule{time}"), el("gallery search")));
    await session.step(33, "And user double-clicks on \"BddCruxMolecule{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BddCruxMolecule{time}\" link in gallery"))));
    await session.step(34, "Then the \"BddCruxMolecule{time}\" view should be current", () => viewIsCurrent(page, session.text("BddCruxMolecule{time}")));
    await session.step(35, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
    await session.step(36, "Then \"BddCruxMolecule{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddCruxMolecule{time}\" dialog")), "visible"));
    await session.step(37, "When user clicks on editor of mol input in \"BddCruxMolecule{time}\" dialog", () => clickOn(page, el(session.text("editor of mol input in \"BddCruxMolecule{time}\" dialog"))));
    await session.step(38, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(39, "And the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(40, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(41, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(42, "And user clicks on the \"atom 1\" area of crux sketcher widget", () => clickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(43, "And user clicks on the \"atom 2\" area of crux sketcher widget", () => clickArea(page, "atom 2", el("crux sketcher widget")));
    await session.step(44, "And user clicks on the \"atom 3\" area of crux sketcher widget", () => clickArea(page, "atom 3", el("crux sketcher widget")));
    await session.step(45, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"CCCCC\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CCCCC"));
    await session.step(46, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(47, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(48, "When user clicks on OK button in \"BddCruxMolecule{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddCruxMolecule{time}\" dialog"))));
    await session.step(49, "Then the \"BddCruxMolecule{time}\" dialog should close", () => dialogCloses(page, session.text("BddCruxMolecule{time}")));
    await session.step(50, "And the script results should show \"got\" as \"CCCCC\"", () => scriptResult(page, "got", "CCCCC"));
  });
});
