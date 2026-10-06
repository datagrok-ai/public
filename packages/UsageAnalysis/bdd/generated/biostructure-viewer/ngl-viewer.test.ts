/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/ngl-viewer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.viewer.ngl]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, selectIn, shouldBe, shouldContainText, shouldOffer, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doesNotOfferColumn} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("NGL viewer on a structure table — empty state and settings", () => {
  const session = feature(test, "features/biostructure-viewer/ngl-viewer.feature", import.meta.url);
  test("NGL viewer on a structure table — empty state and settings", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.viewer.ngl"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(22, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(24, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(25, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(26, "And \"Add viewer\" icon in toolbar should be visible", () => shouldBe(page, el("\"Add viewer\" icon in toolbar"), "visible"));
    await session.step(27, "When user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(28, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(29, "When user types \"NGL\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "NGL", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(30, "And user clicks on first \"NGL\" button in \"Add Viewer\" dialog", () => clickOn(page, el("first \"NGL\" button in \"Add Viewer\" dialog")));
    await session.step(31, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(32, "And NGL viewer should be visible", () => shouldBe(page, el("NGL viewer"), "visible"));
    await session.step(33, "And \"Open...\" link in NGL viewer should be visible", () => shouldBe(page, el("\"Open...\" link in NGL viewer"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("NGL added to a table with only a Molecule3D column shows Open... and no structure", async () => {
      await session.step(38, "Then \"Open...\" link in NGL viewer should be visible", () => shouldBe(page, el("\"Open...\" link in NGL viewer"), "visible"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Style > Representation offers the NGL representations and takes ball+stick", async () => {
      await session.step(42, "When user clicks on settings icon of NGL viewer", () => clickOn(page, el("settings icon of NGL viewer")));
      await session.step(43, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(44, "Then \"Representation\" property in context panel should contain text \"cartoon\"", () => shouldContainText(page, el("\"Representation\" property in context panel"), "cartoon"));
      await session.step(45, "When user clicks on \"Representation\" property in context panel", () => clickOn(page, el("\"Representation\" property in context panel")));
      await session.step(46, "Then \"Representation\" property in context panel should offer \"cartoon, backbone, ball+stick, licorice, hyperball, surface\"", () => shouldOffer(page, el("\"Representation\" property in context panel"), "cartoon, backbone, ball+stick, licorice, hyperball, surface"));
      await session.step(47, "When user selects \"ball+stick\" in \"Representation\" property in context panel", () => selectIn(page, "ball+stick", el("\"Representation\" property in context panel")));
      await session.step(48, "Then \"representation\" property of NGL viewer should be \"ball+stick\"", () => propertyShouldBe(page, "representation", el("NGL viewer"), "ball+stick"));
      await session.step(49, "When user selects \"cartoon\" in \"Representation\" property in context panel", () => selectIn(page, "cartoon", el("\"Representation\" property in context panel")));
      await session.step(50, "Then \"representation\" property of NGL viewer should be \"cartoon\"", () => propertyShouldBe(page, "representation", el("NGL viewer"), "cartoon"));
      await session.step(51, "Given \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(52, "Then \"Ligand\" property in context panel should not offer the column \"pdb\"", () => doesNotOfferColumn(page, el("\"Ligand\" property in context panel"), "pdb"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
