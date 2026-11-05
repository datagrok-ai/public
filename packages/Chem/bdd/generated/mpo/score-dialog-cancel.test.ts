/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/mpo/score-dialog-cancel.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {shippedMpoProfile} from '../../bindings/mpo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, contextPanelOpen, dialogCloses, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MPO Score leaves the table as it was when the dialog is cancelled", () => {
  const session = feature(test, "features/mpo/score-dialog-cancel.feature", import.meta.url);
  test("Cancel after picking a profile leaves the table's own columns only", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the shipped MPO profile \"ADME Test\" is on the server", () => shippedMpoProfile(page, "ADME Test"));
    await session.step(16, "And user opens a table \"adme\" with:", () => openTableOf(page, "adme", [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","-4"],["CCO","-7","3","-5"],["c1ccncc1","-5","1","-3"]]), [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","-4"],["CCO","-7","3","-5"],["c1ccncc1","-5","1","-3"]]);
    await session.step(23, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
    await session.step(25, "Then \"MPO Score\" dialog should be visible", () => shouldBe(page, el("\"MPO Score\" dialog"), "visible"));
    await session.step(26, "When user selects \"+ Create New...\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "+ Create New...", el("Profile input in \"MPO Score\" dialog")));
    await session.step(27, "Then context panel should contain the text \"MPO Profile\"", () => shouldContainText(page, el("context panel"), "MPO Profile"));
    await session.step(28, "When user clicks on CANCEL button in \"MPO Score\" dialog", () => clickOn(page, el("CANCEL button in \"MPO Score\" dialog")));
    await session.step(29, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
    await session.step(30, "And the table should have the columns \"smiles, Caco2, Lipophilicity, Solubility\"", () => columnsExactly(page, "smiles, Caco2, Lipophilicity, Solubility"));
  });
});
