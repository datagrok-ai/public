/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/templates.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reaction-export.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {checksWords, clickDialogButton, clickNamedButton, enterName, enterSmarts, namedButtonDisabled, namedButtonEnabled, namedButtonsStill, namedNoBlock, noteNamedButtons, putBackReactions, reactionSaved, reactionTable, slotShowsWord, slotShowsWords} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-guns in reaction templates: their checks, in the cell's Checks and before OK", () => {
  const session = feature(test, "features/footguns/templates.feature", import.meta.url);
  test("The Checks section of a template cell names what the template does wrong [template=[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]C, level=note, first=no map number, second=new]", {tag: ["@footguns", "@batch-1", "@catalogue-TPL-05"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(20, "Given user opens a table \"templates\" with a reaction column \"template\" of units \"smarts\":", () => reactionTable(page, "templates", "template", "smarts", [["template"],["[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]C"]]), [["template"],["[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]C"]]);
    await session.step(23, "When user clicks on the \"cell 1 of template\" area of grid", () => clickArea(page, "cell 1 of template", el("grid")));
    await session.step(24, "Then the context panel's Checks section lists a note mentioning \"no map number\" and \"new\"", () => checksWords(page, "note", "no map number", "new"));
  });
  test("The Checks section of a template cell names what the template does wrong [template=[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}], level=block, first=placeholder, second=${n}]", {tag: ["@footguns", "@batch-1", "@catalogue-TPL-07"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(20, "Given user opens a table \"templates\" with a reaction column \"template\" of units \"smarts\":", () => reactionTable(page, "templates", "template", "smarts", [["template"],["[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]"]]), [["template"],["[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]"]]);
    await session.step(23, "When user clicks on the \"cell 1 of template\" area of grid", () => clickArea(page, "cell 1 of template", el("grid")));
    await session.step(24, "Then the context panel's Checks section lists a block mentioning \"placeholder\" and \"${n}\"", () => checksWords(page, "block", "placeholder", "${n}"));
  });
  test("The reaction editor blocks SAVE, in Crux mode, on a text RDKit reads and Crux does not", {tag: ["@footguns", "@batch-1", "@catalogue-TPL-17"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(40, "Given user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(43, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(44, "When user picks \"Chem | Transform | Reactions | Transformation...\" from the top menu", () => pickFromTopMenu(page, "Chem | Transform | Reactions | Transformation..."));
    await session.step(45, "And user clicks the \"New Reaction\" button of the dialog", () => clickDialogButton(page, "New Reaction"));
    await session.step(46, "And user enters \"Fragment group test\" as the reaction name", () => enterName(page, "Fragment group test"));
    await session.step(47, "And user enters \"CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|\" as the reaction SMARTS", () => enterSmarts(page, "CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|"));
    await session.step(48, "Then the middle slot of the footer shows a block flag mentioning \"f:\"", () => slotShowsWord(page, "block", "f:"));
    await session.step(49, "And the SAVE button of the \"New Reaction\" dialog is disabled", () => namedButtonDisabled(page, "SAVE", "New Reaction"));
    await session.step(50, "When user enters \"CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:2.3|\" as the reaction SMARTS", () => enterSmarts(page, "CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:2.3|"));
    await session.step(51, "Then the footer of the \"New Reaction\" dialog shows no block flag", () => namedNoBlock(page, "New Reaction"));
    await session.step(52, "And the SAVE button of the \"New Reaction\" dialog is enabled", () => namedButtonEnabled(page, "SAVE", "New Reaction"));
  });
  test("The reaction editor warns and still saves a template whose every run would fail", {tag: ["@footguns", "@batch-1", "@catalogue-TPL-03", "@catalogue-UI-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(58, "Given user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(61, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(62, "And the user's saved reactions are put back when the feature ends", () => putBackReactions(page));
    await session.step(63, "When user picks \"Chem | Transform | Reactions | Transformation...\" from the top menu", () => pickFromTopMenu(page, "Chem | Transform | Reactions | Transformation..."));
    await session.step(64, "And user clicks the \"New Reaction\" button of the dialog", () => clickDialogButton(page, "New Reaction"));
    await session.step(65, "And user notes where the SAVE and CANCEL buttons of the \"New Reaction\" dialog are", () => noteNamedButtons(page, "SAVE", "CANCEL", "New Reaction"));
    await session.step(66, "And user enters \"Duplicate map test\" as the reaction name", () => enterName(page, "Duplicate map test"));
    await session.step(67, "And user enters \"[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]\" as the reaction SMARTS", () => enterSmarts(page, "[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"));
    await session.step(68, "Then the middle slot of the footer shows a warn flag mentioning \"map number 1\" and \"twice\"", () => slotShowsWords(page, "warn", "map number 1", "twice"));
    await session.step(69, "And the SAVE and CANCEL buttons of the \"New Reaction\" dialog are where they were", () => namedButtonsStill(page, "SAVE", "CANCEL", "New Reaction"));
    await session.step(70, "And the SAVE button of the \"New Reaction\" dialog is enabled", () => namedButtonEnabled(page, "SAVE", "New Reaction"));
    await session.step(71, "When user clicks on SAVE button of the \"New Reaction\" dialog", () => clickNamedButton(page, "SAVE", "New Reaction"));
    await session.step(72, "Then the reaction \"Duplicate map test\" is in the user's saved reactions", () => reactionSaved(page, "Duplicate map test"));
  });
});
