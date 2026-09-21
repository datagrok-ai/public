/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/renderers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.rendering, bio.rendering.biln, bio.rendering.separator, bio.rendering.monomer, GROK-12164]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColors, areaPainted, noBalloons, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence cell renderers", () => {
  const session = feature(test, "features/render/renderers.feature", import.meta.url);
  test("Sequence cell renderers", {tag: ["@journey", "@realizes:bio.rendering", "@realizes:bio.rendering.biln", "@realizes:bio.rendering.separator", "@realizes:bio.rendering.monomer", "@realizes:GROK-12164"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("A HELM column renders with the HELM renderer in monomer colors", async () => {
      await session.step(14, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(15, "Then \"HELM string\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM string", "Macromolecule"));
      await session.step(16, "And \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
      await session.step(17, "And the \"cell type of HELM string\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM string", el("grid"), "helm"));
      await session.step(18, "And the \"cell 1 of HELM string\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of HELM string", el("grid"), 3));
      await session.step(19, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A separator column renders with the sequence renderer in monomer colors", async () => {
      await session.step(22, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
      await session.step(23, "Then \"MSA\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "MSA", "Macromolecule"));
      await session.step(24, "And \"MSA\" column should have units \"separator\"", () => columnUnits(page, "MSA", "separator"));
      await session.step(25, "And \"MSA\" column should have tag \"separator\" equal to \"/\"", () => columnTag(page, "MSA", "separator", "/"));
      await session.step(26, "And the \"cell type of MSA\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of MSA", el("grid"), "sequence"));
      await session.step(27, "And the \"cell 1 of MSA\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of MSA", el("grid"), 3));
      await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Converting HELM to separator gives the new column the separator renderer", async () => {
      await session.step(31, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(32, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
      await session.step(33, "And user selects \"separator\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "separator", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
      await session.step(34, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
      await session.step(35, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(36, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(37, "And a new column matching \"^separator\\(HELM string\\)\" should have been added", () => newColumnMatching(page, "^separator\\(HELM string\\)"));
      await session.step(38, "And \"separator(HELM string)\" column should have units \"separator\"", () => columnUnits(page, "separator(HELM string)", "separator"));
      await session.step(39, "And the \"cell type of separator(HELM string)\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of separator(HELM string)", el("grid"), "sequence"));
      await session.step(40, "And \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
      await session.step(41, "And the \"cell type of HELM string\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM string", el("grid"), "helm"));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Split to Monomers columns render with the monomer renderer", async () => {
      await session.step(45, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
      await session.step(46, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
      await session.step(47, "And user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
      await session.step(48, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(49, "And 17 new columns should have been added", () => newColumnsCount(page, 17));
      await session.step(50, "And \"1\" column should have semantic type \"Monomer\"", () => columnSemType(page, "1", "Monomer"));
      await session.step(51, "And the \"cell type of 1\" reading of grid should be \"Monomer\"", () => readingReads(page, "cell type of 1", el("grid"), "Monomer"));
      await session.step(52, "And the \"cell 1 of 1\" area of grid should be painted", () => areaPainted(page, "cell 1 of 1", el("grid")));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
