/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-molecule-parameter.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, currentViewType, packageInstalled, queryTextOnServer, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule, rowMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A query's molecule parameter, sketched in Crux", () => {
  const session = feature(test, "features/queries/query-molecule-parameter.feature", import.meta.url);
  test("A pattern drawn in Crux for the query's Molecule parameter reaches the query as SMILES", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And a query \"BddCruxPattern\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, "BddCruxPattern", "System:Datagrok", "--input: string pattern {semType: Molecule}\nselect @pattern as pattern"));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(27, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(29, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BddCruxPattern tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---Datagrok---BddCruxPattern tree node inside browse tree")));
    await session.step(30, "Then \"BddCruxPattern\" dialog should be visible", () => shouldBe(page, el("\"BddCruxPattern\" dialog"), "visible"));
    await session.step(31, "When user clicks on editor of pattern input in \"BddCruxPattern\" dialog", () => clickOn(page, el("editor of pattern input in \"BddCruxPattern\" dialog")));
    await session.step(32, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(33, "And the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(34, "When user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(35, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(36, "And user clicks on crux nitrogen tool", () => clickOn(page, el("crux nitrogen tool")));
    await session.step(37, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(38, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"c1ccncc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "c1ccncc1"));
    await session.step(39, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(40, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(41, "When user clicks on OK button in \"BddCruxPattern\" dialog", () => clickOn(page, el("OK button in \"BddCruxPattern\" dialog")));
    await session.step(42, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(43, "And the table should have 1 row", () => rowCount(page, 1));
    await session.step(44, "And every value of \"pattern\" column should match \"^\\S+$\"", () => everyValueMatches(page, "pattern", "^\\S+$"));
    await session.step(45, "And the molecule in row 1 of \"pattern\" column should be \"c1ccncc1\"", () => rowMolecule(page, 1, "pattern", "c1ccncc1"));
  });
});
