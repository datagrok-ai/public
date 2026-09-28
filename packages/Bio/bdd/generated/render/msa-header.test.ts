/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/msa-header.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.render.msa-header]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDatasetRowsAs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noBalloons, noErrors, pickFromAreaContextMenu, readingReads, wheelOverAreaTimesHolding} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The WebLogo header of long sequence columns", () => {
  const session = feature(test, "features/render/msa-header.feature", import.meta.url);
  test("The WebLogo header of long sequence columns", {tag: ["@journey", "@realizes:bio.render.msa-header"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(18, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await session.step(20, "And \"AntibodyLC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyLC", "Macromolecule"));
    await run.scenario("Both chains get the header, ruled from 1", async () => {
      await session.step(23, "Then the \"header tracks of AntibodyHC\" reading of grid should be \"Conservation, WebLogo\"", () => readingReads(page, "header tracks of AntibodyHC", el("grid"), "Conservation, WebLogo"));
      await session.step(24, "And the \"header tracks of AntibodyLC\" reading of grid should be \"Conservation, WebLogo\"", () => readingReads(page, "header tracks of AntibodyLC", el("grid"), "Conservation, WebLogo"));
      await session.step(25, "And the \"header positions of AntibodyHC\" reading of grid should be \"1, 10, 20\"", () => readingReads(page, "header positions of AntibodyHC", el("grid"), "1, 10, 20"));
      await session.step(26, "And the \"header positions of AntibodyLC\" reading of grid should be \"1, 10, 20\"", () => readingReads(page, "header positions of AntibodyLC", el("grid"), "1, 10, 20"));
    });
    await run.scenario("Extract Region cuts CDR3 out of the aligned column with the scheme's positions", async () => {
      await session.step(29, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(30, "And user clicks on OK button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("OK button in \"Apply Antibody Numbering\" dialog")));
      await session.step(31, "Then a new column \"AntibodyHC (aligned)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)"));
      await session.step(32, "And the \"header tracks of AntibodyHC (aligned)\" reading of grid should be \"Conservation, WebLogo, Annotations\"", () => readingReads(page, "header tracks of AntibodyHC (aligned)", el("grid"), "Conservation, WebLogo, Annotations"));
      await session.step(33, "When user picks \"Bio > Calculate > Extract Region...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Extract Region..."));
      await session.step(34, "And user selects \"AntibodyHC (aligned)\" in Sequence input in \"Get Sequence Region\" dialog", () => selectIn(page, "AntibodyHC (aligned)", el("Sequence input in \"Get Sequence Region\" dialog")));
      await session.step(35, "And user selects \"CDR3: 105-117\" in Region input in \"Get Sequence Region\" dialog", () => selectIn(page, "CDR3: 105-117", el("Region input in \"Get Sequence Region\" dialog")));
      await session.step(36, "And user enters \"CDR3\" into \"Column name\" input in \"Get Sequence Region\" dialog", () => enterInto(page, "CDR3", el("\"Column name\" input in \"Get Sequence Region\" dialog")));
      await session.step(37, "And user clicks on OK button in \"Get Sequence Region\" dialog", () => clickOn(page, el("OK button in \"Get Sequence Region\" dialog")));
      await session.step(38, "Then a new column \"CDR3\" should have been added", () => newColumnNamed(page, "CDR3"));
      await session.step(39, "And \"CDR3\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "CDR3", "aligned", "SEQ.MSA"));
      await session.step(40, "When user scrolls the mouse wheel down 5 times over the \"row header 1\" area of grid holding Shift", () => wheelOverAreaTimesHolding(page, "down", 5, "row header 1", el("grid"), "Shift"));
      await session.step(41, "Then the \"header positions of CDR3\" reading of grid should be \"105, 110\"", () => readingReads(page, "header positions of CDR3", el("grid"), "105, 110"));
      await session.step(42, "When user hovers over the \"position 110 of CDR3 header\" area of grid", () => hoverArea(page, "position 110 of CDR3 header", el("grid")));
      await session.step(43, "Then tooltip should contain text \"Position: 110\"", () => shouldContainText(page, el("tooltip"), "Position: 110"));
      await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The aligned cell's Extract CDR3 keeps the scheme's positions too", async () => {
      await session.step(47, "When user scrolls the mouse wheel up 5 times over the \"row header 1\" area of grid holding Shift", () => wheelOverAreaTimesHolding(page, "up", 5, "row header 1", el("grid"), "Shift"));
      await session.step(48, "And user picks \"Annotations > Extract CDR3 as Column\" from the context menu of the \"cell 1 of AntibodyHC (aligned)\" area of grid", () => pickFromAreaContextMenu(page, "Annotations > Extract CDR3 as Column", "cell 1 of AntibodyHC (aligned)", el("grid")));
      await session.step(49, "Then a new column \"AntibodyHC (aligned)(CDR3)\" should have been added", () => newColumnNamed(page, "AntibodyHC (aligned)(CDR3)"));
      await session.step(50, "When user scrolls the mouse wheel down 5 times over the \"row header 1\" area of grid holding Shift", () => wheelOverAreaTimesHolding(page, "down", 5, "row header 1", el("grid"), "Shift"));
      await session.step(51, "Then the \"header positions of AntibodyHC (aligned)(CDR3)\" reading of grid should be \"105, 110\"", () => readingReads(page, "header positions of AntibodyHC (aligned)(CDR3)", el("grid"), "105, 110"));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
