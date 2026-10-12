/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/hit-design-crux.feature
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
import {campaignOnServer, campaignSaved} from '../../bindings/hit-design.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openApp, packageInstalled, sketcherIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule, rowMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A row added in Hit Design, drawn in Crux", () => {
  const session = feature(test, "features/apps/hit-design-crux.feature", import.meta.url);
  test("A new row opens Crux empty and ready, and a molecule drawn there is written to the row with its V-iD", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(14, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Hit Design campaign \"BDD Crux\" is on the server, and is as it was again when the feature ends", () => campaignOnServer(page, "BDD Crux"));
    await session.step(20, "Given user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await session.step(21, "When user clicks on \"BDD Crux\" link", () => clickOn(page, el("\"BDD Crux\" link")));
    await session.step(22, "Then the \"BDD Crux\" view should be current", () => viewIsCurrent(page, "BDD Crux"));
    await session.step(23, "And the table should have 1 row", () => rowCount(page, 1));
    await session.step(24, "When user clicks on \"Add new row\" icon in toolbar", () => clickOn(page, el("\"Add new row\" icon in toolbar")));
    await session.step(25, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(26, "And the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(27, "And the \"atoms\" reading of crux sketcher widget should be 0", () => readingIs(page, "atoms", el("crux sketcher widget"), 0));
    await session.step(28, "When user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(29, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(30, "And user clicks on crux nitrogen tool", () => clickOn(page, el("crux nitrogen tool")));
    await session.step(31, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(32, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"c1ccncc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "c1ccncc1"));
    await session.step(33, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(34, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(35, "And the table should have 2 rows", () => rowCount(page, 2));
    await session.step(36, "And the molecule in row 2 of \"Molecule\" column should be \"c1ccncc1\"", () => rowMolecule(page, 2, "Molecule", "c1ccncc1"));
    await session.step(37, "And the Hit Design campaign \"BDD Crux\" should be saved with 2 rows, row 2 the molecule \"c1ccncc1\" with its V-iD", () => campaignSaved(page, "BDD Crux", 2, 2, "c1ccncc1"));
  });
});
