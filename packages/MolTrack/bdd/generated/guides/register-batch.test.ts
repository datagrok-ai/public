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
import {clickOn, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register a batch of a registered compound", () => {
  const session = feature(test, "features/guides/register-batch.feature", import.meta.url);
  test("Open the compound, add a batch to it, save it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(12, "And no batch with the EPA Batch ID \"GUIDE-001\" is registered, now and at feature end", () => noBatch(page, "GUIDE-001"));
    await session.step(14, "And user opens the \"Compounds?corporate_compound_id=DG-000482\" view of MolTrack EMS", () => openEmsView(page, "Compounds?corporate_compound_id=DG-000482"), undefined, "Open the compound by its Grok ID: DG-000482, bisphenol A");
    enter(page, "MolTrack EMS");
    await session.step(15, "Then breadcrumbs should contain text \"DG-000482\"", () => shouldContainText(page, el("breadcrumbs"), "DG-000482"));
    await session.step(17, "When user clicks on \"New\" button in batches tab", () => clickOn(page, el("\"New\" button in batches tab")), undefined, "Its batches are below it; New starts one with the compound filled in");
    await session.step(19, "And user types \"GUIDE-001\" into \"EPA Batch ID\" text input", () => typeInto(page, "GUIDE-001", el("\"EPA Batch ID\" text input")), undefined, "The batch's own fields");
    await session.step(20, "And user types \"98\" into \"Purity (%)\" number input", () => typeInto(page, "98", el("\"Purity (%)\" number input")));
    await session.step(21, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(23, "Then 1 batch with the EPA Batch ID \"GUIDE-001\" should be on the server", () => batchCount(page, 1, "GUIDE-001"), undefined, "The batch is saved on that compound");
  });
});
