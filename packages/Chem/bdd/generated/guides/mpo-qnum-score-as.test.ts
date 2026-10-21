/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/mpo-qnum-score-as.feature
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
import {clickOn, selectIn, shouldBe, shouldHaveValue, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, dialogCloses, openTableOf, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Choose how MPO scores values with < or >", () => {
  const session = feature(test, "features/guides/mpo-qnum-score-as.feature", import.meta.url);
  test("Score the bounds of a qnum property by their limit value", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the shipped MPO profile \"ADME Test\" is on the server", () => shippedMpoProfile(page, "ADME Test"));
    await session.step(21, "And user opens a table \"ADME assays\" with:", () => openTableOf(page, "ADME assays", [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","<-5","2","0"],["CCO",">-5","2","0"],["c1ccncc1","-6","2","0"],["CC(=O)O","<-9","2","0"]]), [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","<-5","2","0"],["CCO",">-5","2","0"],["c1ccncc1","-6","2","0"],["CC(=O)O","<-9","2","0"]]);
    await session.step(27, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
    await session.step(28, "And user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
    await session.step(29, "And user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
    await session.step(30, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
    await session.step(31, "And a new column \"MPO ADME Test\" should have been added", () => newColumnNamed(page, "MPO ADME Test"));
    await session.step(32, "And the \"MPO ADME Test\" cell of row 1 should be displayed as \"0.71429\"", () => displayedInRow(page, "MPO ADME Test", 1, "0.71429"));
    await session.step(33, "And the \"MPO ADME Test\" cell of row 2 should be displayed as \"0.71429\"", () => displayedInRow(page, "MPO ADME Test", 2, "0.71429"));
    await session.step(34, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
    await session.step(35, "And user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
    await session.step(36, "And user switches on \"Design mode\" input in \"MPO Score\" dialog", () => switchOn(page, el("\"Design mode\" input in \"MPO Score\" dialog")));
    await session.step(37, "And user clicks on \"Settings\" icon in first MPO property", () => clickOn(page, el("\"Settings\" icon in first MPO property")));
    await session.step(38, "Then \"Caco2\" dialog should be visible", () => shouldBe(page, el("\"Caco2\" dialog"), "visible"));
    await session.step(39, "And \"Score as\" input in \"Caco2\" dialog should have value \"Least favorable value\"", () => shouldHaveValue(page, el("\"Score as\" input in \"Caco2\" dialog"), "Least favorable value"));
    await session.step(40, "When user selects \"Limit value\" in \"Score as\" input in \"Caco2\" dialog", () => selectIn(page, "Limit value", el("\"Score as\" input in \"Caco2\" dialog")));
    await session.step(41, "And user clicks on OK button in \"Caco2\" dialog", () => clickOn(page, el("OK button in \"Caco2\" dialog")));
    await session.step(42, "Then the \"Caco2\" dialog should close", () => dialogCloses(page, "Caco2"));
    await session.step(43, "When user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
    await session.step(44, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
    await session.step(45, "And the \"MPO ADME Test\" cell of row 1 should be displayed as \"0.85714\"", () => displayedInRow(page, "MPO ADME Test", 1, "0.85714"));
    await session.step(46, "And the \"MPO ADME Test\" cell of row 2 should be displayed as \"0.85714\"", () => displayedInRow(page, "MPO ADME Test", 2, "0.85714"));
    await session.step(47, "And the \"MPO ADME Test\" cell of row 3 should be displayed as \"1.00000\"", () => displayedInRow(page, "MPO ADME Test", 3, "1.00000"));
    await session.step(48, "And the \"MPO ADME Test\" cell of row 4 should be displayed as \"0.71429\"", () => displayedInRow(page, "MPO ADME Test", 4, "0.71429"));
  });
});
