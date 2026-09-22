/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/molecule-rendering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.molecule-rendering-end-to-end, chem.int.render-feeds-search, GROK-16870]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnCount, columnSemType, everyValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColors, areaNotColor, areaPainted, hoverArea, noBalloons, noErrors, oneTooltip, painted, tooltipSomeColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Molecule, reaction and mixture cells are drawn, in the grid and in a viewer tooltip", () => {
  const session = feature(test, "features/render/molecule-rendering.feature", import.meta.url);
  test("Molecule, reaction and mixture cells are drawn, in the grid and in a viewer tooltip", {tag: ["@journey", "@realizes:chem.cp.molecule-rendering-end-to-end", "@realizes:chem.int.render-feeds-search", "@realizes:GROK-16870"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("Molecule cells are drawn as structures, not as text", async () => {
      await session.step(21, "Then \"canonical_smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "canonical_smiles", "Molecule"));
      await session.step(22, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(23, "And the \"cell 1 of canonical_smiles\" area of grid should be painted", () => areaPainted(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(24, "And the \"cell 1 of canonical_smiles\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 1 of canonical_smiles", el("grid"), 1));
      await session.step(25, "And the \"cell 2 of canonical_smiles\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 2 of canonical_smiles", el("grid"), 1));
      await session.step(26, "And the \"cell 3 of canonical_smiles\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 3 of canonical_smiles", el("grid"), 1));
      await session.step(27, "And the \"cell 1 of molregno\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "cell 1 of molregno", el("grid"), "#FF0000"));
      await session.step(28, "And the \"cell 1 of molregno\" area of grid should not contain the color \"#0000FF\"", () => areaNotColor(page, "cell 1 of molregno", el("grid"), "#0000FF"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A scatter plot tooltip draws the molecule of the row under the pointer", async () => {
      await session.step(32, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["x","NumHAcceptors"],["y","NumHDonors"]]), [["x","NumHAcceptors"],["y","NumHDonors"]]);
      await session.step(35, "Then scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(36, "When user hovers over the \"marker of row 1\" area of scatter plot viewer", () => hoverArea(page, "marker of row 1", el("scatter plot viewer")));
      await session.step(37, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(38, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(39, "And the canvases of tooltip should be painted in at least 2 colors", () => canvasColors(page, el("tooltip"), 2));
      await session.step(40, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A box plot tooltip draws a molecule and the viewer keeps painting", async () => {
      await session.step(45, "Given user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["value","NumHAcceptors"]]), [["value","NumHAcceptors"]]);
      await session.step(47, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(48, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(49, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(50, "And the canvases of tooltip should be painted in at least 2 colors", () => canvasColors(page, el("tooltip"), 2));
      await session.step(51, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Reaction cells are typed and drawn by the reaction renderer", async () => {
      await session.step(56, "Given user opens test-reactions dataset", () => openDataset(page, ds("test-reactions")));
      await session.step(57, "Then the table should have 17 rows", () => rowCount(page, 17));
      await session.step(58, "And \"reaction\" column should have semantic type \"ChemicalReaction\"", () => columnSemType(page, "reaction", "ChemicalReaction"));
      await session.step(59, "And every value of \"reaction\" column should contain \">>\"", () => everyValueContains(page, "reaction", ">>"));
      await session.step(60, "And the \"cell 1 of reaction\" area of grid should be painted", () => areaPainted(page, "cell 1 of reaction", el("grid")));
      await session.step(61, "And the \"cell 1 of reaction\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 1 of reaction", el("grid"), 1));
      await session.step(62, "And the \"cell 2 of reaction\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 2 of reaction", el("grid"), 1));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Mixture cells are typed and drawn by the mixture renderer", async () => {
      await session.step(67, "Given user opens test_mixtures dataset", () => openDataset(page, ds("test_mixtures")));
      await session.step(68, "Then the table should have 4 rows", () => rowCount(page, 4));
      await session.step(69, "And the table should have 1 column", () => columnCount(page, 1));
      await session.step(70, "And \"mixture\" column should have semantic type \"ChemicalMixture\"", () => columnSemType(page, "mixture", "ChemicalMixture"));
      await session.step(71, "And every value of \"mixture\" column should contain \"mixfileVersion\"", () => everyValueContains(page, "mixture", "mixfileVersion"));
      await session.step(72, "And every value of \"mixture\" column should contain \"contents\"", () => everyValueContains(page, "mixture", "contents"));
      await session.step(73, "And the \"cell 1 of mixture\" area of grid should be painted", () => areaPainted(page, "cell 1 of mixture", el("grid")));
      await session.step(74, "And the \"cell 1 of mixture\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 1 of mixture", el("grid"), 1));
      await session.step(75, "And the \"cell 2 of mixture\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 2 of mixture", el("grid"), 1));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
