/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/18b-composition.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {gridScrolledToRow} from '../../bindings/demo-search.js';
import {bioInitialized, selectedByMonomer} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hoverArea, noBalloons, noErrors, painted, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Composition", () => {
  const session = feature(test, "features/demos/18b-composition.feature", import.meta.url);
  test("Pick out the peptides with one residue at one position", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(23, "And the \"sequence\" column of the grid is 300 px wide", () => gridColumnWidth(page, "sequence", 300));
    await session.step(24, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","activity"]]), [["valueColumnName","activity"]]);
    await session.step(26, "And the grid is scrolled down to row 52", () => gridScrolledToRow(page, 52));
    await session.step(27, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
    await session.step(28, "Then WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(29, "Given user sets \"positionHeight\" property of WebLogo viewer to \"100%\"", () => setProperty(page, "positionHeight", el("WebLogo viewer"), "100%"));
    await session.step(30, "Then WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
    await session.step(31, "When the video says \"Click a residue at a position: its peptides light up everywhere\"", () => videoSays(page, "Click a residue at a position: its peptides light up everywhere"));
    await session.step(32, "And user hovers over the \"monomer N at position 3\" area of WebLogo viewer", () => hoverArea(page, "monomer N at position 3", el("WebLogo viewer")));
    await session.step(33, "And user clicks on the \"monomer N at position 3\" area of WebLogo viewer", () => clickArea(page, "monomer N at position 3", el("WebLogo viewer")));
    await session.step(34, "Then 31 rows should be selected", () => selectedRowCount(page, 31));
    await session.step(35, "And only the rows with \"N\" at position 3 of \"sequence\" column should be selected", () => selectedByMonomer(page, "N", 3, "sequence"));
    await session.step(36, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
    await session.step(37, "When user hovers over the \"monomer V at position 14\" area of WebLogo viewer", () => hoverArea(page, "monomer V at position 14", el("WebLogo viewer")));
    await session.step(38, "And the video says \"Another residue, another family: here the least active\"", () => videoSays(page, "Another residue, another family: here the least active"));
    await session.step(39, "And user clicks on the \"monomer V at position 14\" area of WebLogo viewer", () => clickArea(page, "monomer V at position 14", el("WebLogo viewer")));
    await session.step(40, "Then 34 rows should be selected", () => selectedRowCount(page, 34));
    await session.step(41, "And only the rows with \"V\" at position 14 of \"sequence\" column should be selected", () => selectedByMonomer(page, "V", 14, "sequence"));
    await session.step(42, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
