/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/reaction-cell-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.reactions.cell-editor]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-export.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cellAsRemembered, cellDiffers, cruxKey, rememberCell} from '../../bindings/crux.js';
import {engineGate} from '../../bindings/engine.js';
import {editorArranged, editorArrows, readsAsReaction, scrollGridToRow} from '../../bindings/reaction-editor.js';
import {openProjectCsv} from '../../bindings/reactions.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {doubleClickArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A reaction cell edited in the Crux sketcher", () => {
  const session = feature(test, "features/sketcher/reaction-cell-editor.feature", import.meta.url);
  test("A reaction cell edited in the Crux sketcher", {tag: ["@journey", "@realizes:chem.reactions.cell-editor"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(12, "And user opens the \"fixtures/BDD-reactions-showcase.csv\" file of the project", () => openProjectCsv(page, "fixtures/BDD-reactions-showcase.csv"));
    await run.scenario("A double-click opens the cell's reaction drawn; OK without an edit leaves it byte for byte", async () => {
      await session.step(17, "When user remembers the value of \"reaction\" column in row 1", () => rememberCell(page, "reaction", 1));
      await session.step(18, "And user double-clicks on the \"cell 1 of reaction\" area of grid", () => doubleClickArea(page, "cell 1 of reaction", el("grid")));
      await session.step(19, "Then reaction editor dialog should be visible", () => shouldBe(page, el("reaction editor dialog"), "visible"));
      await session.step(20, "And the reaction editor should draw 1 reaction arrow", () => editorArrows(page, 1));
      await session.step(21, "When user clicks on OK button in reaction editor dialog", () => clickOn(page, el("OK button in reaction editor dialog")));
      await session.step(22, "Then reaction editor dialog should be absent", () => shouldBe(page, el("reaction editor dialog"), "absent"));
      await session.step(23, "And the value of \"reaction\" column in row 1 should be, byte for byte, as remembered", () => cellAsRemembered(page, "reaction", 1));
      await session.step(24, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An edit, confirmed with OK, is written into the cell, still a reaction", async () => {
      await session.step(28, "When user remembers the value of \"reaction\" column in row 1", () => rememberCell(page, "reaction", 1));
      await session.step(29, "And user picks \"Edit reaction...\" from the context menu of the \"cell 1 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 1 of reaction", el("grid")));
      await session.step(30, "Then reaction editor dialog should be visible", () => shouldBe(page, el("reaction editor dialog"), "visible"));
      await session.step(31, "When user presses the \"Delete\" key over the \"atom 0\" area of Crux sketcher widget", () => cruxKey(page, "Delete", "atom 0", el("Crux sketcher widget")));
      await session.step(32, "And user clicks on OK button in reaction editor dialog", () => clickOn(page, el("OK button in reaction editor dialog")));
      await session.step(33, "Then reaction editor dialog should be absent", () => shouldBe(page, el("reaction editor dialog"), "absent"));
      await session.step(34, "And the value of \"reaction\" column in row 1 should differ from what was remembered", () => cellDiffers(page, "reaction", 1));
      await session.step(35, "And Chem should read the value of \"reaction\" column in row 1 as a reaction", () => readsAsReaction(page, "reaction", 1));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A route cell opens arranged in rows; OK without an edit leaves it byte for byte", async () => {
      await session.step(40, "When user remembers the value of \"reaction\" column in row 10", () => rememberCell(page, "reaction", 10));
      await session.step(41, "And user scrolls the grid to row 10 of the \"reaction\" column", () => scrollGridToRow(page, 10, "reaction"));
      await session.step(42, "And user picks \"Edit reaction...\" from the context menu of the \"cell 10 of reaction\" area of grid", () => pickFromAreaContextMenu(page, "Edit reaction...", "cell 10 of reaction", el("grid")));
      await session.step(43, "Then reaction editor dialog should be visible", () => shouldBe(page, el("reaction editor dialog"), "visible"));
      await session.step(44, "And the reaction editor should draw the reaction in at least 2 rows, its bonds at least 16 px long", () => editorArranged(page, 2, 16));
      await session.step(45, "When user clicks on OK button in reaction editor dialog", () => clickOn(page, el("OK button in reaction editor dialog")));
      await session.step(46, "Then reaction editor dialog should be absent", () => shouldBe(page, el("reaction editor dialog"), "absent"));
      await session.step(47, "And the value of \"reaction\" column in row 10 should be, byte for byte, as remembered", () => cellAsRemembered(page, "reaction", 10));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
