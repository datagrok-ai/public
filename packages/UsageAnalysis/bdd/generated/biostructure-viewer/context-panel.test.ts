/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/context-panel.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.panel.3d-structure, biostructureviewer.panel.pdb-information]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {collapse, doubleClickOn, isExpanded, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, contextPanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Context panel of a Molecule3D cell — 3D Structure and PDB Information", () => {
  const session = feature(test, "features/biostructure-viewer/context-panel.feature", import.meta.url);
  test("Context panel of a Molecule3D cell — 3D Structure and PDB Information", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.panel.3d-structure", "@realizes:biostructureviewer.panel.pdb-information"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(25, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(26, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(27, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(28, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(29, "Given the context panel is open", () => contextPanelOpen(page));
    await run.scenario("3D Structure embeds a Biostructure viewer for the current Molecule3D cell", async () => {
      await session.step(32, "When user clicks on the \"cell 1 of pdb\" area of grid", () => clickArea(page, "cell 1 of pdb", el("grid")));
      await session.step(33, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(34, "Given \"3D Structure\" section in context panel is expanded", () => isExpanded(page, el("\"3D Structure\" section in context panel")));
      await session.step(35, "Then \"Reset Camera\" button in \"3D Structure\" section in context panel should be visible", () => shouldBe(page, el("\"Reset Camera\" button in \"3D Structure\" section in context panel"), "visible"));
      await session.step(36, "Given \"PDB Information\" section in context panel is expanded", () => isExpanded(page, el("\"PDB Information\" section in context panel")));
      await session.step(37, "Then \"PDB Information\" section in context panel should contain text \"ASPARTYL PROTEASE\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "ASPARTYL PROTEASE"));
      await session.step(38, "When user clicks on the \"cell 2 of pdb\" area of grid", () => clickArea(page, "cell 2 of pdb", el("grid")));
      await session.step(39, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(40, "And \"PDB Information\" section in context panel should contain text \"HYDROLASE\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "HYDROLASE"));
      await session.step(41, "And \"PDB Information\" section in context panel should not contain text \"ASPARTYL PROTEASE\"", () => shouldNotContainText(page, el("\"PDB Information\" section in context panel"), "ASPARTYL PROTEASE"));
      await session.step(42, "And \"Reset Camera\" button in \"3D Structure\" section in context panel should be visible", () => shouldBe(page, el("\"Reset Camera\" button in \"3D Structure\" section in context panel"), "visible"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
      await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(45, "When user collapses \"3D Structure\" section in context panel", () => collapse(page, el("\"3D Structure\" section in context panel")));
      await session.step(46, "Then no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("PDB Information reads the header of the current cell's PDB text", async () => {
      await session.step(50, "When user clicks on the \"cell 1 of pdb\" area of grid", () => clickArea(page, "cell 1 of pdb", el("grid")));
      await session.step(51, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(52, "And \"PDB Information\" section in context panel should contain text \"ASPARTYL PROTEASE\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "ASPARTYL PROTEASE"));
      await session.step(53, "And \"PDB Information\" section in context panel should contain text \"HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY"));
      await session.step(54, "And \"PDB Information\" section in context panel should contain text \"rcsb.org/structure/1QBS\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "rcsb.org/structure/1QBS"));
      await session.step(55, "When user clicks on the \"cell 2 of pdb\" area of grid", () => clickArea(page, "cell 2 of pdb", el("grid")));
      await session.step(56, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(57, "And \"PDB Information\" section in context panel should contain text \"HYDROLASE\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "HYDROLASE"));
      await session.step(58, "And \"PDB Information\" section in context panel should contain text \"HIV PROTEASE WITH INHIBITOR AB-2\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "HIV PROTEASE WITH INHIBITOR AB-2"));
      await session.step(59, "And \"PDB Information\" section in context panel should contain text \"rcsb.org/structure/1ZP8\"", () => shouldContainText(page, el("\"PDB Information\" section in context panel"), "rcsb.org/structure/1ZP8"));
      await session.step(60, "And \"PDB Information\" section in context panel should not contain text \"ASPARTYL PROTEASE\"", () => shouldNotContainText(page, el("\"PDB Information\" section in context panel"), "ASPARTYL PROTEASE"));
      await session.step(61, "When user collapses \"PDB Information\" section in context panel", () => collapse(page, el("\"PDB Information\" section in context panel")));
      await session.step(62, "Then no errors should have been logged", () => noErrors(page));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
