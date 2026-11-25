/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/reaction-cell-editor-lossy.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.reactions.cell-editor]
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
import {namedButtonEnabled, reactionTable, slotShowsWord, valueContains} from '../../bindings/footguns.js';
import {valueLacks, writeOverArrow} from '../../bindings/reactions-sketch-stack.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A reaction cell edited with text over its arrow warns before OK writes it, and never writes \"*\"", () => {
  const session = feature(test, "features/sketcher/reaction-cell-editor-lossy.feature", import.meta.url);
  test("A reaction cell edited with text over its arrow warns before OK writes it, and never writes \"*\"", {tag: ["@journey", "@realizes:chem.reactions.cell-editor", "@eval-b7"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await run.scenario("Text over the arrow is named in the footer before OK, and OK writes the reaction without a \"*\"", async () => {
      await session.step(14, "Given user opens a table \"lossy\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "lossy", "reaction", "smiles", [["reaction"],["CC(=O)O.CCN>>CCNC(C)=O"]]), [["reaction"],["CC(=O)O.CCN>>CCNC(C)=O"]]);
      await session.step(17, "When user picks \"Edit reaction...\" from the context menu of the \"cell 1 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 1 of reaction", el("grid")));
      await session.step(18, "Then reaction editor dialog should be visible", () => shouldBe(page, el("reaction editor dialog"), "visible"));
      await session.step(19, "When user writes \"HATU, DIPEA\" over arrow 1 in the reaction editor", () => writeOverArrow(page, "HATU, DIPEA", 1));
      await session.step(20, "Then the middle slot of the footer shows a warn flag mentioning \"texts\"", () => slotShowsWord(page, "warn", "texts"));
      await session.step(21, "And the OK button of the \"Reaction\" dialog is enabled", () => namedButtonEnabled(page, "OK", "Reaction"));
      await session.step(22, "When user clicks on OK button in reaction editor dialog", () => clickOn(page, el("OK button in reaction editor dialog")));
      await session.step(23, "Then reaction editor dialog should be absent", () => shouldBe(page, el("reaction editor dialog"), "absent"));
      await session.step(24, "And the value in row 1 of \"reaction\" column does not contain \"*\"", () => valueLacks(page, 1, "reaction", "*"));
      await session.step(25, "And the value in row 1 of \"reaction\" column contains \">>\"", () => valueContains(page, 1, "reaction", ">>"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
