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
import {clickOn, pressKeyIn, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register a batch of a registered compound", () => {
  const session = feature(test, "features/guides/register-batch.feature", import.meta.url);
  test("Pick the compound, fill in the batch, save it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(12, "And no batch with the EPA Batch ID \"GUIDE-001\" is registered, now and at feature end", () => noBatch(page, "GUIDE-001"));
    await session.step(13, "And user opens the \"Batches\" view of MolTrack EMS", () => openEmsView(page, "Batches"));
    enter(page, "MolTrack EMS");
    await session.step(15, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Start a new batch");
    await session.step(17, "And user types \"Bisphenol A\" into \"Compound\" domain pick", () => typeInto(page, "Bisphenol A", el("\"Compound\" domain pick")), undefined, "Pick its compound by Grok ID or name");
    await session.step(18, "And user presses Enter in \"Compound\" domain pick", () => pressKeyIn(page, "Enter", el("\"Compound\" domain pick")));
    await session.step(20, "Then compound card should contain text \"Bisphenol A\"", () => shouldContainText(page, el("compound card"), "Bisphenol A"), undefined, "The compound the batch is of");
    await session.step(22, "When user types \"GUIDE-001\" into \"EPA Batch ID\" text input", () => typeInto(page, "GUIDE-001", el("\"EPA Batch ID\" text input")), undefined, "The batch's own fields");
    await session.step(23, "And user types \"98.2\" into \"Purity\" number input", () => typeInto(page, "98.2", el("\"Purity\" number input")));
    await session.step(24, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(26, "Then notification should contain text \"saved\"", () => shouldContainText(page, el("notification"), "saved"), undefined, "The batch is saved on that compound");
    await session.step(27, "And compound card should contain text \"Bisphenol A\"", () => shouldContainText(page, el("compound card"), "Bisphenol A"));
    await session.step(28, "And 1 batch with the EPA Batch ID \"GUIDE-001\" should be on the server", () => batchCount(page, 1, "GUIDE-001"));
  });
});
