/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/reaction-editor.feature
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
import {clickReactionEditorButton, corpusRxnTable, eraseArrow, methylOnFirstReactant, methylOnProductMap, noReactantFlag, noteReactionEditorButtons, reactionEditorButtonDisabled, reactionEditorButtonEnabled, reactionEditorButtonsStill, reactionEditorSlotEmpty, reactionTable, rxnNotationWarn, slotShowsWord, slotShowsWords, undoEdit, valueAsOpened, valueStillContains} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-guns in the reaction cell editor (a reaction cell's Edit reaction...)", () => {
  const session = feature(test, "features/footguns/reaction-editor.feature", import.meta.url);
  test("OK and Cancel of the reaction cell editor stay put when a flag comes and goes; erasing the arrow blocks OK", {tag: ["@footguns", "@batch-3", "@needs-land-5", "@catalogue-UI-01", "@catalogue-RXN-16"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(24, "Given user opens a table \"reactions\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "reactions", "reaction", "smiles", [["reaction"],["CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"]]), [["reaction"],["CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"]]);
    await session.step(27, "When user picks \"Edit reaction...\" from the context menu of the \"cell 1 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 1 of reaction", el("grid")));
    await session.step(28, "Then the footer of the reaction editor has an empty middle slot", () => reactionEditorSlotEmpty(page));
    await session.step(29, "And user notes where the OK and CANCEL buttons of the reaction editor are", () => noteReactionEditorButtons(page, "OK", "CANCEL"));
    await session.step(30, "When user erases the arrow in the reaction editor", () => eraseArrow(page));
    await session.step(31, "Then the middle slot of the footer shows a block flag mentioning \"arrow\"", () => slotShowsWord(page, "block", "arrow"));
    await session.step(32, "And the OK button of the reaction editor is disabled", () => reactionEditorButtonDisabled(page, "OK"));
    await session.step(33, "And the OK and CANCEL buttons of the reaction editor are where they were", () => reactionEditorButtonsStill(page, "OK", "CANCEL"));
    await session.step(34, "When user undoes the last edit in the reaction editor", () => undoEdit(page));
    await session.step(35, "Then the footer of the reaction editor has an empty middle slot", () => reactionEditorSlotEmpty(page));
    await session.step(36, "And the OK and CANCEL buttons of the reaction editor are where they were", () => reactionEditorButtonsStill(page, "OK", "CANCEL"));
    await session.step(37, "When user clicks on CANCEL button of the reaction editor", () => clickReactionEditorButton(page, "CANCEL"));
    await session.step(38, "Then the value in row 1 of \"reaction\" column is exactly what it was", () => valueAsOpened(page, 1, "reaction"));
  });
  test("A value the editor reads as a molecule is not written back as one [units=smiles, value=CC >> CO, keyword=space]", {tag: ["@footguns", "@batch-3", "@needs-land-5", "@catalogue-RXN-15"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(42, "Given user opens a table \"reactions\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "reactions", "reaction", "smiles", [["reaction"],["CC >> CO"]]), [["reaction"],["CC >> CO"]]);
    await session.step(45, "When user picks \"Edit reaction...\" from the context menu of the \"cell 1 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 1 of reaction", el("grid")));
    await session.step(46, "Then the middle slot of the footer shows a block flag mentioning \"reaction\" and \"space\"", () => slotShowsWords(page, "block", "reaction", "space"));
    await session.step(47, "And the OK button of the reaction editor is disabled", () => reactionEditorButtonDisabled(page, "OK"));
    await session.step(48, "When user clicks on CANCEL button of the reaction editor", () => clickReactionEditorButton(page, "CANCEL"));
    await session.step(49, "Then the value in row 1 of \"reaction\" column is exactly what it was", () => valueAsOpened(page, 1, "reaction"));
  });
  test("Editing one component of a template keeps the others as written and names any query feature lost", {tag: ["@footguns", "@batch-3", "@needs-land-5", "@catalogue-TPL-12"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(57, "Given user opens a table \"templates\" with a reaction column \"template\" of units \"smarts\":", () => reactionTable(page, "templates", "template", "smarts", [["template"],["[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]"]]), [["template"],["[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]"]]);
    await session.step(60, "When user picks \"Edit reaction...\" from the context menu of the \"cell 1 of template\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 1 of template", el("grid")));
    await session.step(61, "And user adds a methyl to the product's carbon 2 in the reaction editor", () => methylOnProductMap(page, 2));
    await session.step(62, "Then the footer shows no flag about the reactants", () => noReactantFlag(page));
    await session.step(63, "When user clicks on OK button of the reaction editor", () => clickReactionEditorButton(page, "OK"));
    await session.step(64, "Then the value in row 1 of \"template\" column still contains \"[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1]\"", () => valueStillContains(page, 1, "template", "[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1]"));
    await session.step(65, "And the value in row 1 of \"template\" column still contains \"[OH,O-]\"", () => valueStillContains(page, 1, "template", "[OH,O-]"));
  });
  test("An RXN block edited in the editor warns that it is written in another notation", {tag: ["@footguns", "@batch-3", "@needs-land-5", "@catalogue-RXN-17"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(69, "Given user opens a table \"rxn\" with a reaction column \"reaction\" of units \"rxn\" holding the corpus's \"block-v2000\" RXN block", () => corpusRxnTable(page, "rxn", "reaction", "rxn", "block-v2000"));
    await session.step(70, "When user picks \"Edit reaction...\" from the context menu of the \"cell 1 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 1 of reaction", el("grid")));
    await session.step(71, "And user adds a methyl to the first reactant in the reaction editor", () => methylOnFirstReactant(page));
    await session.step(72, "Then the middle slot of the footer shows a warn flag mentioning \"RXN\" and the notation it writes instead", () => rxnNotationWarn(page, "RXN"));
    await session.step(73, "And the OK button of the reaction editor is enabled", () => reactionEditorButtonEnabled(page, "OK"));
  });
});
