/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/register-batch.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {batchCount, demoData, noBatch, openEmsView} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldContainText, typeAtCaret} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack: register a batch of a registered compound", () => {
  const session = feature(test, "features/guides/register-batch.feature", import.meta.url);
  test("Open the compound, add a batch to it, save it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the MolTrack demo data is on the server", () => demoData(page));
    await session.step(12, "And no batch with the Source Batch Code \"GUIDE-001\" is registered, now and at feature end", () => noBatch(page, "GUIDE-001"));
    await session.step(14, "And user opens the \"Compounds?corporate_compound_id=DG-000482\" view of MolTrack", () => openEmsView(page, "Compounds?corporate_compound_id=DG-000482"), undefined, "Open the compound by its Grok ID: DG-000482, bisphenol A");
    enter(page, "MolTrack");
    await session.step(15, "Then breadcrumbs should contain text \"DG-000482\"", () => shouldContainText(page, el("breadcrumbs"), "DG-000482"));
    await session.step(17, "When user clicks on \"New\" button in batches tab", () => clickOn(page, el("\"New\" button in batches tab")), undefined, "Its batches are below it; New starts one with the compound filled in");
    await session.step(19, "And user clicks on the \"cell 2 of source_batch_code\" area of grid in batches tab", () => clickArea(page, "cell 2 of source_batch_code", el("grid in batches tab")), undefined, "The batch's own fields, typed into its new row under bisphenol A's one batch");
    await session.step(20, "And user types \"GUIDE-001\" at the caret", () => typeAtCaret(page, "GUIDE-001"));
    await session.step(21, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(22, "And user clicks on the \"cell 2 of purity\" area of grid in batches tab", () => clickArea(page, "cell 2 of purity", el("grid in batches tab")));
    await session.step(23, "And user types \"98\" at the caret", () => typeAtCaret(page, "98"));
    await session.step(24, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(25, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(27, "Then 1 batch with the Source Batch Code \"GUIDE-001\" should be on the server", () => batchCount(page, 1, "GUIDE-001"), undefined, "The batch is saved on that compound");
  });
});
