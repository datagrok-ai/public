/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/weblogo-residue-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized, selectedByMonomer} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, openToolbox, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pick out the peptides with one residue at one position", () => {
  const session = feature(test, "features/guides/weblogo-residue-selection.feature", import.meta.url);
  test("Click a letter of the WebLogo to select the peptides that have it", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
    await session.step(18, "Then WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(19, "When user opens toolbox", () => openToolbox(page));
    await session.step(20, "And user clicks on histogram icon on toolbox", () => clickOn(page, el("histogram icon on toolbox")));
    await session.step(21, "And user selects \"activity\" in Value column input in histogram viewer", () => selectIn(page, "activity", el("Value column input in histogram viewer")));
    await session.step(22, "And user clicks on the \"monomer N at position 3\" area of WebLogo viewer", () => clickArea(page, "monomer N at position 3", el("WebLogo viewer")));
    await session.step(23, "Then 31 rows should be selected", () => selectedRowCount(page, 31));
    await session.step(24, "And only the rows with \"N\" at position 3 of \"sequence\" column should be selected", () => selectedByMonomer(page, "N", 3, "sequence"));
    await session.step(25, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
  });
});
