/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/flags-footer.feature
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
import {buttonsClearOfSlot, checksNone, checksWords, hoverSlot, noteSketcherButtons, putMolecule, reactionTable, resizeSketcherDialog, sketcherButtonsStill, slotCut, slotEmptySketcher, slotShows, slotTooltip} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, doubleClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-gun flags: where they show, and that OK and Cancel never move", () => {
  const session = feature(test, "features/footguns/flags-footer.feature", import.meta.url);
  test("A flag longer than the slot is cut and shows whole on hover, and the buttons do not move", {tag: ["@footguns", "@batch-1", "@catalogue-UI-02", "@catalogue-UI-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(24, "Given user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(27, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(28, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(29, "And user notes where the OK and CANCEL buttons of the sketcher dialog are", () => noteSketcherButtons(page, "OK", "CANCEL"));
    await session.step(30, "And user puts \"C(C)(C)(C)(C)C.N(C)(C)(C)(C)C.c1cccc1.[CH3:1]C(C)(C)(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "C(C)(C)(C)(C)C.N(C)(C)(C)(C)C.c1cccc1.[CH3:1]C(C)(C)(C)(C)C"));
    await session.step(31, "Then the middle slot of the footer shows a flag whose text is cut to the slot's width", () => slotCut(page));
    await session.step(32, "And the OK and CANCEL buttons of the sketcher dialog are where they were", () => sketcherButtonsStill(page, "OK", "CANCEL"));
    await session.step(33, "When user hovers over the middle slot of the footer", () => hoverSlot(page));
    await session.step(34, "Then a tooltip shows the flag's whole text, naming every problem found", () => slotTooltip(page));
    await session.step(35, "When user puts \"CCN\" into the Crux sketcher of the dialog", () => putMolecule(page, "CCN"));
    await session.step(36, "Then the footer of the sketcher dialog has an empty middle slot", () => slotEmptySketcher(page));
    await session.step(37, "And the OK and CANCEL buttons of the sketcher dialog are where they were", () => sketcherButtonsStill(page, "OK", "CANCEL"));
  });
  test("A dialog narrowed by the user keeps its buttons in place when a flag appears", {tag: ["@footguns", "@batch-1", "@catalogue-UI-02"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(41, "Given user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(44, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(45, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(46, "And user resizes the sketcher dialog to 420 px wide", () => resizeSketcherDialog(page, 420));
    await session.step(47, "And user notes where the OK and CANCEL buttons of the sketcher dialog are", () => noteSketcherButtons(page, "OK", "CANCEL"));
    await session.step(48, "And user puts \"C(C)(C)(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "C(C)(C)(C)(C)C"));
    await session.step(49, "Then the middle slot of the footer shows a block flag", () => slotShows(page, "block"));
    await session.step(50, "And the OK and CANCEL buttons of the sketcher dialog are where they were", () => sketcherButtonsStill(page, "OK", "CANCEL"));
    await session.step(51, "And the OK and CANCEL buttons do not overlap the middle slot", () => buttonsClearOfSlot(page));
  });
  test("The Checks section lists the current cell's flags and follows the current cell", {tag: ["@footguns", "@batch-1", "@catalogue-UI-04"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(55, "Given user opens a table \"templates\" with a reaction column \"template\" of units \"smarts\":", () => reactionTable(page, "templates", "template", "smarts", [["template"],["[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]"]]), [["template"],["[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]"]]);
    await session.step(59, "When user clicks on the \"cell 1 of template\" area of grid", () => clickArea(page, "cell 1 of template", el("grid")));
    await session.step(60, "Then the context panel's Checks section lists a warn mentioning \"map number 1\" and \"twice\"", () => checksWords(page, "warn", "map number 1", "twice"));
    await session.step(61, "When user clicks on the \"cell 2 of template\" area of grid", () => clickArea(page, "cell 2 of template", el("grid")));
    await session.step(62, "Then the context panel's Checks section lists no flag", () => checksNone(page));
    await session.step(63, "When user clicks on the \"cell 1 of template\" area of grid", () => clickArea(page, "cell 1 of template", el("grid")));
    await session.step(64, "Then the context panel's Checks section lists a warn mentioning \"map number 1\" and \"twice\"", () => checksWords(page, "warn", "map number 1", "twice"));
  });
});
