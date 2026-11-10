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
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterToAnyOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · BILN connections", () => {
  const session = feature(test, "features/demos/02-biln-connections.feature", import.meta.url);
  test("Rings, cross-links and chains drawn as connections", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "And user opens BILN_complex dataset", () => openDataset(page, ds("BILN_complex")));
    await session.step(18, "Then \"biln\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "biln", "Macromolecule"));
    await session.step(19, "And \"biln\" column should have units \"biln\"", () => columnUnits(page, "biln", "biln"));
    await session.step(20, "Given user filters rows where \"name\" is one of \"oxytocin, linaclotide, kalata_B1, cyclosporin_like, stapled_helix, insulin_human, MAP4_NANP, trivalent_cRGDfK_THA\"", () => filterToAnyOf(page, "name", "oxytocin, linaclotide, kalata_B1, cyclosporin_like, stapled_helix, insulin_human, MAP4_NANP, trivalent_cRGDfK_THA"));
    await session.step(21, "And the grid shows only the columns \"name, category, biln\"", () => gridShowsColumns(page, "name, category, biln"));
    await session.step(22, "And the \"biln\" column of the grid is 1040 px wide", () => gridColumnWidth(page, "biln", 1040));
    await session.step(23, "And the monomers of the \"biln\" column are drawn at 16 px", () => sequenceFontSize(page, "biln", 16));
    await session.step(24, "And the rows of the grid are 86 px tall", () => gridRowHeight(page, 86));
    await session.step(25, "When user hovers over monomer 1 of the \"cell 1 of biln\" area of grid", () => hoverMonomer(page, 1, "cell 1 of biln", el("grid")));
    await session.step(26, "And the video says \"Disulfides, rings and cross-links: brackets over the monomers\"", () => videoSays(page, "Disulfides, rings and cross-links: brackets over the monomers"));
    await session.step(27, "And user hovers over monomer 5 of the \"cell 4 of biln\" area of grid", () => hoverMonomer(page, 5, "cell 4 of biln", el("grid")));
    await session.step(28, "And user hovers over monomer 1 of the \"cell 33 of biln\" area of grid", () => hoverMonomer(page, 1, "cell 33 of biln", el("grid")));
    await session.step(29, "And user hovers over monomer 4 of the \"cell 42 of biln\" area of grid", () => hoverMonomer(page, 4, "cell 42 of biln", el("grid")));
    await session.step(30, "And the video says \"Two chains, one line each, bonded across\"", () => videoSays(page, "Two chains, one line each, bonded across"));
    await session.step(31, "And user hovers over monomer 7 of the \"cell 48 of biln\" area of grid", () => hoverMonomer(page, 7, "cell 48 of biln", el("grid")));
    await session.step(32, "And user hovers over monomer 9 of the \"cell 61 of biln\" area of grid", () => hoverMonomer(page, 9, "cell 61 of biln", el("grid")));
    await session.step(33, "And user hovers over monomer 16 of the \"cell 62 of biln\" area of grid", () => hoverMonomer(page, 16, "cell 62 of biln", el("grid")));
    await session.step(34, "Then no errors should have been logged", () => noErrors(page));
  });
});
