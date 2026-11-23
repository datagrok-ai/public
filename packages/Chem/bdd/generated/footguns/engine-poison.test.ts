/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/engine-poison.feature
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
import {engineStops, enterInSketcherDialog, moleculeTable, okDisabledOrNot, slotShowsWord, valueAsOpened} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A Crux sketcher whose engine stops, or that cannot write its drawing, never writes over the cell", () => {
  const session = feature(test, "features/footguns/engine-poison.feature", import.meta.url);
  test("The engine stopping under the molecule cell editor warns quietly, and neither Enter nor OK empties the cell", {tag: ["@footguns", "@engine-poison"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(20, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"CCO\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "CCO"));
    await session.step(21, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(22, "And the Crux engine of the sketcher dialog stops", () => engineStops(page));
    await session.step(23, "Then the middle slot of the footer shows a warn flag mentioning \"could not be read back\"", () => slotShowsWord(page, "warn", "could not be read back"));
    await session.step(24, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(25, "When user presses Enter while the sketcher dialog has the focus", () => enterInSketcherDialog(page));
    await session.step(26, "Then the value in row 1 of \"molecule\" column is exactly what it was", () => valueAsOpened(page, 1, "molecule"));
    await session.step(27, "When user clicks the OK button of the sketcher dialog, disabled or not", () => okDisabledOrNot(page));
    await session.step(28, "Then the value in row 1 of \"molecule\" column is exactly what it was", () => valueAsOpened(page, 1, "molecule"));
  });
});
