/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/biostructure-viewer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.biostructure]
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
import {clickOn, doubleClickOn, downloadContains, fileDownloaded, isExpanded, selectIn, shouldBe, shouldContainText, typeInto, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu, propertyShouldBe, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Biostructure viewer on a structure table — structure column, representation, camera, download", () => {
  const session = feature(test, "features/biostructure-viewer/biostructure-viewer.feature", import.meta.url);
  test("Biostructure viewer on a structure table — structure column, representation, camera, download", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.viewer.biostructure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(32, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(33, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(34, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(35, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(36, "And \"pdb\" column should have semantic type \"Molecule3D\"", () => columnSemType(page, "pdb", "Molecule3D"));
    await session.step(37, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
    await session.step(38, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(39, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(40, "When user types \"Biostructure\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "Biostructure", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(41, "And user clicks on first \"Biostructure\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Biostructure\" button in \"Add Viewer\" dialog")));
    await session.step(42, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(43, "And Biostructure viewer should be visible", () => shouldBe(page, el("Biostructure viewer"), "visible"));
    await session.step(44, "When user clicks on settings icon of Biostructure viewer", () => clickOn(page, el("settings icon of Biostructure viewer")));
    await session.step(45, "And user selects \"pdb\" in \"Biostructure Id\" property in context panel", () => selectIn(page, "pdb", el("\"Biostructure Id\" property in context panel")));
    await session.step(46, "Then \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
    await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("The viewer shows the pdb column's structure with the cartoon representation", async () => {
      await session.step(51, "Then \"Data File\" input in Biostructure viewer should be absent", () => shouldBe(page, el("\"Data File\" input in Biostructure viewer"), "absent"));
      await session.step(52, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(53, "Then \"Representation\" property in context panel should contain text \"cartoon\"", () => shouldContainText(page, el("\"Representation\" property in context panel"), "cartoon"));
      await session.step(54, "When user clicks on the \"cell 3 of id\" area of grid", () => clickArea(page, "cell 3 of id", el("grid")));
      await session.step(55, "Then row 3 should be current", () => currentRowIs(page, 3));
      await session.step(56, "And \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The representation switches from the settings (GROK-11759)", async () => {
      await session.step(61, "When user selects \"ball-and-stick\" in \"Representation\" property in context panel", () => selectIn(page, "ball-and-stick", el("\"Representation\" property in context panel")));
      await session.step(62, "Then \"representation\" property of Biostructure viewer should be \"ball-and-stick\"", () => propertyShouldBe(page, "representation", el("Biostructure viewer"), "ball-and-stick"));
      await session.step(63, "And the \"representation\" reading of Biostructure viewer should be \"ball-and-stick\"", () => readingReads(page, "representation", el("Biostructure viewer"), "ball-and-stick"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(66, "When user selects \"molecular-surface\" in \"Representation\" property in context panel", () => selectIn(page, "molecular-surface", el("\"Representation\" property in context panel")));
      await session.step(67, "Then \"representation\" property of Biostructure viewer should be \"molecular-surface\"", () => propertyShouldBe(page, "representation", el("Biostructure viewer"), "molecular-surface"));
      await session.step(68, "And the \"representation\" reading of Biostructure viewer should be \"molecular-surface\"", () => readingReads(page, "representation", el("Biostructure viewer"), "molecular-surface"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(71, "When user selects \"cartoon\" in \"Representation\" property in context panel", () => selectIn(page, "cartoon", el("\"Representation\" property in context panel")));
      await session.step(72, "Then \"representation\" property of Biostructure viewer should be \"cartoon\"", () => propertyShouldBe(page, "representation", el("Biostructure viewer"), "cartoon"));
      await session.step(73, "And the \"representation\" reading of Biostructure viewer should be \"cartoon\"", () => readingReads(page, "representation", el("Biostructure viewer"), "cartoon"));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
      await session.step(75, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Reset Camera overlay button keeps the viewer working", async () => {
      await session.step(78, "When user clicks on \"Reset Camera\" button in Biostructure viewer", () => clickOn(page, el("\"Reset Camera\" button in Biostructure viewer")));
      await session.step(79, "Then \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
      await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The viewer's context menu downloads the structure as PDB and as CIF", async () => {
      await session.step(84, "Given user watches downloads", () => watchDownloads(page));
      await session.step(85, "When user opens the context menu of Biostructure viewer", () => openContextMenu(page, el("Biostructure viewer")));
      await session.step(86, "Then the open menu should list \"Download > As PDB\"", () => menuLists(page, "Download > As PDB"));
      await session.step(87, "And the open menu should list \"Download > As CIF\"", () => menuLists(page, "Download > As CIF"));
      await session.step(88, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(89, "And user picks \"Download > As PDB\" from the context menu of Biostructure viewer", () => pickFromContextMenu(page, "Download > As PDB", el("Biostructure viewer")));
      await session.step(90, "Then a file \"file.pdb\" should have been downloaded", () => fileDownloaded(page, "file.pdb"));
      await session.step(91, "And the downloaded file \"file.pdb\" should contain text \"CRYST1   42.100   54.600   69.000\"", () => downloadContains(page, "file.pdb", "CRYST1   42.100   54.600   69.000"));
      await session.step(92, "When user picks \"Download > As CIF\" from the context menu of Biostructure viewer", () => pickFromContextMenu(page, "Download > As CIF", el("Biostructure viewer")));
      await session.step(93, "Then a file \"file.cif\" should have been downloaded", () => fileDownloaded(page, "file.cif"));
      await session.step(94, "And the downloaded file \"file.cif\" should contain text \"_atom_site\"", () => downloadContains(page, "file.cif", "_atom_site"));
      await session.step(95, "And the downloaded file \"file.cif\" should contain text \"PROTO-ONCOGENE TYROSINE-PROTEIN KINASE SRC\"", () => downloadContains(page, "file.cif", "PROTO-ONCOGENE TYROSINE-PROTEIN KINASE SRC"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
      await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
