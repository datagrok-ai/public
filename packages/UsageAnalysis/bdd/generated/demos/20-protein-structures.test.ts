/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/20-protein-structures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelClosed, browsePanelOpen, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridRowHeight, gridShowsColumns, viewersSettled} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("6 · Antibodies and proteins · Protein structures", () => {
  const session = feature(test, "features/demos/20-protein-structures.feature", import.meta.url);
  test("Protein structures in 3D", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(22, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(24, "When the video says \"Structure files preview in 3D, right from Browse\"", () => videoSays(page, "Structure files preview in 3D, right from Browse"));
    await session.step(25, "And user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree")));
    await session.step(26, "Then the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(27, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(28, "When user clicks on Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree")));
    await session.step(29, "Then the \"1RQ9.mmcif\" view should be current", () => viewIsCurrent(page, "1RQ9.mmcif"));
    await session.step(30, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(31, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(32, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(33, "And \"pdb\" column should have semantic type \"Molecule3D\"", () => columnSemType(page, "pdb", "Molecule3D"));
    await session.step(34, "Given the browse panel is closed", () => browsePanelClosed(page));
    await session.step(35, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(36, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(37, "And the grid shows only the columns \"pdb_id, prop2, desc\"", () => gridShowsColumns(page, "pdb_id, prop2, desc"));
    await session.step(38, "And the rows of the grid are 32 px tall", () => gridRowHeight(page, 32));
    await session.step(39, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(40, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(41, "When user types \"Biostructure\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Biostructure", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(42, "And user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
    await session.step(43, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(44, "And Biostructure viewer should be visible", () => shouldBe(page, el("Biostructure viewer"), "visible"));
    await session.step(45, "Given the viewers have settled", () => viewersSettled(page));
    await session.step(46, "When user clicks on the \"cell 3 of pdb_id\" area of grid", () => clickArea(page, "cell 3 of pdb_id", el("grid")));
    await session.step(47, "Then row 3 should be current", () => currentRowIs(page, 3));
    await session.step(48, "Given the viewers have settled", () => viewersSettled(page));
    await session.step(49, "When the video says \"The viewer follows the current row\"", () => videoSays(page, "The viewer follows the current row"));
    await session.step(50, "And user clicks on the \"cell 6 of pdb_id\" area of grid", () => clickArea(page, "cell 6 of pdb_id", el("grid")));
    await session.step(51, "Then row 6 should be current", () => currentRowIs(page, 6));
    await session.step(52, "Given the viewers have settled", () => viewersSettled(page));
    await session.step(53, "Then no errors should have been logged", () => noErrors(page));
  });
});
