/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/03-notation-conversion.feature
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
import {clickOn, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Notation conversion", () => {
  const session = feature(test, "features/demos/03-notation-conversion.feature", import.meta.url);
  test("HELM to BILN and back to plain text", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(14, "And user opens HELM_cyclic dataset", () => openDataset(page, ds("HELM_cyclic")));
    await session.step(15, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(16, "And \"sequence\" column should have units \"helm\"", () => columnUnits(page, "sequence", "helm"));
    await session.step(17, "Given the grid shows only the columns \"sequence_id, sequence\"", () => gridShowsColumns(page, "sequence_id, sequence"));
    await session.step(18, "And the \"sequence\" column of the grid is 300 px wide", () => gridColumnWidth(page, "sequence", 300));
    await session.step(19, "And the rows of the grid are 130 px tall", () => gridRowHeight(page, 130));
    await session.step(20, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
    await session.step(21, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
    await session.step(22, "When user selects \"biln\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "biln", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
    await session.step(23, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(24, "Then a new column \"biln(sequence)\" should have been added", () => newColumnNamed(page, "biln(sequence)"));
    await session.step(25, "And \"biln(sequence)\" column should have units \"biln\"", () => columnUnits(page, "biln(sequence)", "biln"));
    await session.step(26, "Given the \"biln(sequence)\" column of the grid is 560 px wide", () => gridColumnWidth(page, "biln(sequence)", 560));
    await session.step(27, "And the monomers of the \"biln(sequence)\" column are drawn at 14 px", () => sequenceFontSize(page, "biln(sequence)", 14));
    await session.step(28, "When the video says \"In BILN the head-to-tail bond is part of the line, drawn as a bracket\"", () => videoSays(page, "In BILN the head-to-tail bond is part of the line, drawn as a bracket"));
    await session.step(29, "And user hovers over monomer 1 of the \"cell 1 of biln(sequence)\" area of grid", () => hoverMonomer(page, 1, "cell 1 of biln(sequence)", el("grid")));
    await session.step(30, "And user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
    await session.step(31, "And user selects \"separator\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "separator", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
    await session.step(32, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(33, "Then a new column \"separator(sequence)\" should have been added", () => newColumnNamed(page, "separator(sequence)"));
    await session.step(34, "Given the \"separator(sequence)\" column of the grid is 560 px wide", () => gridColumnWidth(page, "separator(sequence)", 560));
    await session.step(35, "And the monomers of the \"separator(sequence)\" column are drawn at 14 px", () => sequenceFontSize(page, "separator(sequence)", 14));
    await session.step(36, "Then no errors should have been logged", () => noErrors(page));
  });
});
