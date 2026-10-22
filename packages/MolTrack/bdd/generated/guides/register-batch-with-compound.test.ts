/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/register-batch-with-compound.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {batchCount, compoundCount, noCompound, openEmsView} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register a batch together with its compound", () => {
  const session = feature(test, "features/guides/register-batch-with-compound.feature", import.meta.url);
  test("A new batch registers its new compound in the same save, a second batch joins it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And no compound with the structure \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" is registered, now and at feature end", () => noCompound(page, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"));
    await session.step(12, "And user opens the \"Batches\" view of MolTrack EMS", () => openEmsView(page, "Batches"));
    enter(page, "MolTrack EMS");
    await session.step(14, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "Start a new batch: the form for its compound comes with it");
    await session.step(15, "Then new compound section should be visible", () => shouldBe(page, el("new compound section"), "visible"));
    await session.step(16, "When user clicks on molecule canvas in new compound section", () => clickOn(page, el("molecule canvas in new compound section")));
    await session.step(17, "And user types \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" into molecule input of sketcher dialog", () => typeInto(page, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO", el("molecule input of sketcher dialog")));
    await session.step(18, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(19, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(20, "And user types \"Octacosane-1,28-diol\" into \"Common Name\" text input in new compound section", () => typeInto(page, "Octacosane-1,28-diol", el("\"Common Name\" text input in new compound section")));
    await session.step(22, "And user types \"GUIDE-001\" into \"EPA Batch ID\" text input", () => typeInto(page, "GUIDE-001", el("\"EPA Batch ID\" text input")), undefined, "The batch's own fields");
    await session.step(23, "And user types \"98.2\" into \"Purity\" number input", () => typeInto(page, "98.2", el("\"Purity\" number input")));
    await session.step(24, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(26, "Then notification should contain text \"saved\"", () => shouldContainText(page, el("notification"), "saved"), undefined, "One save registered the compound and its batch");
    await session.step(27, "And 1 compound with the structure \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" should be registered on the server", () => compoundCount(page, 1, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"));
    await session.step(28, "And the compound with the structure \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" should have 1 batch on the server", () => batchCount(page, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO", 1));
    await session.step(30, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")), undefined, "A second batch of the same structure");
    await session.step(31, "And user clicks on molecule canvas in new compound section", () => clickOn(page, el("molecule canvas in new compound section")));
    await session.step(32, "And user types \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" into molecule input of sketcher dialog", () => typeInto(page, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO", el("molecule input of sketcher dialog")));
    await session.step(33, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(34, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(36, "Then registered compound note should contain text \"the batch is added to it\"", () => shouldContainText(page, el("registered compound note"), "the batch is added to it"), undefined, "The structure is registered: the batch is added to that compound");
    await session.step(37, "And new compound section should be hidden", () => shouldBe(page, el("new compound section"), "hidden"));
    await session.step(38, "When user types \"GUIDE-002\" into \"EPA Batch ID\" text input", () => typeInto(page, "GUIDE-002", el("\"EPA Batch ID\" text input")));
    await session.step(39, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
    await session.step(40, "Then 1 compound with the structure \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" should be registered on the server", () => compoundCount(page, 1, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"));
    await session.step(41, "And the compound with the structure \"OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO\" should have 2 batches on the server", () => batchCount(page, "OCCCCCCCCCCCCCCCCCCCCCCCCCCCCO", 2));
  });
});
