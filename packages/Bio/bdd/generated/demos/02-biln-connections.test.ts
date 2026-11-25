/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/02-biln-connections.feature
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
import {monomerNameLength} from '../../bindings/demo-notation.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterPasses, filterToAnyOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("1 · Any notation · BILN connections", () => {
  const session = feature(test, "features/demos/02-biln-connections.feature", import.meta.url);
  test("Rings, cross-links and conjugates drawn as connections", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(25, "And user opens BILN_complex dataset", () => openDataset(page, ds("BILN_complex")));
    await session.step(26, "Then \"biln\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "biln", "Macromolecule"));
    await session.step(27, "And \"biln\" column should have units \"biln\"", () => columnUnits(page, "biln", "biln"));
    await session.step(28, "Given the browse panel is closed", () => browsePanelClosed(page));
    await session.step(29, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(30, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(31, "And user filters rows where \"name\" is one of \"oxytocin, linaclotide, kalata_B1, cyclosporin_like, stapled_helix, insulin_human, MAP4_NANP, trivalent_cRGDfK_THA\"", () => filterToAnyOf(page, "name", "oxytocin, linaclotide, kalata_B1, cyclosporin_like, stapled_helix, insulin_human, MAP4_NANP, trivalent_cRGDfK_THA"));
    await session.step(32, "And the grid shows only the columns \"name, category, biln\"", () => gridShowsColumns(page, "name, category, biln"));
    await session.step(33, "And the \"biln\" column of the grid is 1040 px wide", () => gridColumnWidth(page, "biln", 1040));
    await session.step(34, "And the monomers of the \"biln\" column are drawn at 16 px", () => sequenceFontSize(page, "biln", 16));
    await session.step(35, "And the monomer names of the \"biln\" column are cut at 8 characters", () => monomerNameLength(page, "biln", 8));
    await session.step(36, "And the rows of the grid are 86 px tall", () => gridRowHeight(page, 86));
    await session.step(37, "When the video says \"Disulfides, rings and cross-links: brackets over the monomers\"", () => videoSays(page, "Disulfides, rings and cross-links: brackets over the monomers"));
    await session.step(38, "And user hovers over monomer 6 of the \"cell 1 of biln\" area of grid", () => hoverMonomer(page, 6, "cell 1 of biln", el("grid")));
    await session.step(39, "And user hovers over monomer 29 of the \"cell 33 of biln\" area of grid", () => hoverMonomer(page, 29, "cell 33 of biln", el("grid")));
    await session.step(40, "And user hovers over monomer 11 of the \"cell 35 of biln\" area of grid", () => hoverMonomer(page, 11, "cell 35 of biln", el("grid")));
    await session.step(41, "And user hovers over monomer 14 of the \"cell 42 of biln\" area of grid", () => hoverMonomer(page, 14, "cell 42 of biln", el("grid")));
    await session.step(42, "And the video says \"Two chains, one line each, bonded across\"", () => videoSays(page, "Two chains, one line each, bonded across"));
    await session.step(43, "And user hovers over monomer 40 of the \"cell 48 of biln\" area of grid", () => hoverMonomer(page, 40, "cell 48 of biln", el("grid")));
    await session.step(44, "And user hovers over monomer 15 of the \"cell 62 of biln\" area of grid", () => hoverMonomer(page, 15, "cell 62 of biln", el("grid")));
    await session.step(45, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(46, "And the video says \"Conjugates: separate chains, bonded side by side\"", () => videoSays(page, "Conjugates: separate chains, bonded side by side"));
    await session.step(47, "Given user filters rows where \"name\" is one of \"Fmoc-FF_hydrogelator, quencher_linker_Lys, R8-bAla-PNA, peptide-drug_conjugate_like\"", () => filterToAnyOf(page, "name", "Fmoc-FF_hydrogelator, quencher_linker_Lys, R8-bAla-PNA, peptide-drug_conjugate_like"));
    await session.step(48, "And the rows of the grid are 96 px tall", () => gridRowHeight(page, 96));
    await session.step(49, "Then 4 rows should pass the filter", () => filterPasses(page, 4));
    await session.step(50, "When user hovers over monomer 2 of the \"cell 64 of biln\" area of grid", () => hoverMonomer(page, 2, "cell 64 of biln", el("grid")));
    await session.step(51, "And user hovers over monomer 6 of the \"cell 68 of biln\" area of grid", () => hoverMonomer(page, 6, "cell 68 of biln", el("grid")));
    await session.step(52, "And user hovers over monomer 9 of the \"cell 71 of biln\" area of grid", () => hoverMonomer(page, 9, "cell 71 of biln", el("grid")));
    await session.step(53, "And user hovers over monomer 17 of the \"cell 75 of biln\" area of grid", () => hoverMonomer(page, 17, "cell 75 of biln", el("grid")));
    await session.step(54, "Then no errors should have been logged", () => noErrors(page));
  });
});
