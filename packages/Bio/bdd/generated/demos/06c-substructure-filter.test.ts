/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/06c-substructure-filter.feature
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
import {noCurrentRow} from '../../bindings/demo-atomic.js';
import {atomicLevelDefaults, gridSortedBy} from '../../bindings/demo-chem.js';
import {monomerNameLength} from '../../bindings/demo-notation.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDatasetRows, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("2 · Every sequence is a molecule · Substructure filter", () => {
  const session = feature(test, "features/demos/06c-substructure-filter.feature", import.meta.url);
  test("One chemical motif across every peptide", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And user opens BILN_complex dataset keeping the first 63 rows", () => openDatasetRows(page, ds("BILN_complex"), 63));
    await session.step(23, "And the \"biln\" column is converted to atomic level", () => atomicLevelDefaults(page, "biln"));
    await session.step(24, "And the grid is sorted by the \"name\" column in descending order", () => gridSortedBy(page, "name", "descending"));
    await session.step(25, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(27, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(28, "And the grid shows only the columns \"name, biln, molfile(biln)\"", () => gridShowsColumns(page, "name, biln, molfile(biln)"));
    await session.step(29, "And the \"biln\" column of the grid is 500 px wide", () => gridColumnWidth(page, "biln", 500));
    await session.step(30, "And the \"molfile(biln)\" column of the grid is 560 px wide", () => gridColumnWidth(page, "molfile(biln)", 560));
    await session.step(31, "And the monomers of the \"biln\" column are drawn at 14 px", () => sequenceFontSize(page, "biln", 14));
    await session.step(32, "And the monomer names of the \"biln\" column are cut at 8 characters", () => monomerNameLength(page, "biln", 8));
    await session.step(33, "And the rows of the grid are 120 px tall", () => gridRowHeight(page, 120));
    await session.step(34, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(35, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(36, "And user adds a card for \"molfile(biln)\" to the filter panel", () => addCardFor(page, "molfile(biln)"));
    await session.step(37, "When user clicks on \"Sketch\" text in \"molfile(biln)\" filter card", () => clickOn(page, el("\"Sketch\" text in \"molfile(biln)\" filter card")));
    await session.step(38, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(39, "When user types \"CSSC\" into molecule input of sketcher dialog", () => typeInto(page, "CSSC", el("molecule input of sketcher dialog")));
    await session.step(40, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(41, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(42, "Then 44 rows should pass the filter", () => filterPasses(page, 44));
    await session.step(43, "When the video says \"44 of 63 peptides carry a disulfide bridge\"", () => videoSays(page, "44 of 63 peptides carry a disulfide bridge"));
    await session.step(44, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
