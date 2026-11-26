/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/03-composition-with-a-weblogo.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioReady} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaAtLeastTall, clickArea, hoverArea, hoverFirstArea, noBalloons, noErrors, painted, someHighlight, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/03-composition-with-a-weblogo.feature", import.meta.url);
  test("See the composition with a WebLogo", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(18, "And the Bio package is initialized", () => bioReady(page));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(22, "And the \"sequence\" column of the grid is 300 px wide", () => gridColumnWidth(page, "sequence", 300));
    await session.step(23, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","activity"]]), [["valueColumnName","activity"]]);
    await session.step(26, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"), undefined, "Bio > Analyze > Composition");
    await session.step(27, "Then WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(28, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
    await session.step(29, "When the video says \"A WebLogo: the residues at each position, the most common on top\" at the top", () => videoSaysAt(page, "A WebLogo: the residues at each position, the most common on top", "top"));
    await session.step(30, "And user hovers over the \"monomer N at position 3\" area of WebLogo viewer", () => hoverArea(page, "monomer N at position 3", el("WebLogo viewer")));
    await session.step(31, "And the video says \"Click a letter to select the peptides that have it there\" at the top", () => videoSaysAt(page, "Click a letter to select the peptides that have it there", "top"));
    await session.step(32, "And user clicks on the \"monomer N at position 3\" area of WebLogo viewer", () => clickArea(page, "monomer N at position 3", el("WebLogo viewer")));
    await session.step(33, "Then 31 rows should be selected", () => selectedRowCount(page, 31));
    await session.step(34, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
    await session.step(36, "When user scrolls the mouse wheel down 3 times over the \"cell 10 of sequence\" area of grid", () => wheelOverAreaTimes(page, "down", 3, "cell 10 of sequence", el("grid")), undefined, "Scroll the table to them");
    await session.step(37, "Then the \"row header 40\" area of grid should be at least 10 pixels tall", () => areaAtLeastTall(page, "row header 40", el("grid"), 10));
    await session.step(38, "When the video says \"Selected everywhere: in the histogram, all near the top of the activity\" at the top", () => videoSaysAt(page, "Selected everywhere: in the histogram, all near the top of the activity", "top"));
    await session.step(39, "And user hovers over the \"monomer V at position 14\" area of WebLogo viewer", () => hoverArea(page, "monomer V at position 14", el("WebLogo viewer")));
    await session.step(40, "And the video says \"Another residue, another family: nearly all at the low end\" at the top", () => videoSaysAt(page, "Another residue, another family: nearly all at the low end", "top"));
    await session.step(41, "And user clicks on the \"monomer V at position 14\" area of WebLogo viewer", () => clickArea(page, "monomer V at position 14", el("WebLogo viewer")));
    await session.step(42, "Then 34 rows should be selected", () => selectedRowCount(page, 34));
    await session.step(43, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
    await session.step(44, "When user scrolls the mouse wheel down 3 times over the \"cell 45 of sequence\" area of grid", () => wheelOverAreaTimes(page, "down", 3, "cell 45 of sequence", el("grid")));
    await session.step(45, "Then the \"row header 70\" area of grid should be at least 10 pixels tall", () => areaAtLeastTall(page, "row header 70", el("grid"), 10));
    await session.step(46, "When user hovers over the first \"selected bin\" area of histogram viewer", () => hoverFirstArea(page, "selected bin", el("histogram viewer")));
    await session.step(47, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
