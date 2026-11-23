/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/leaks/leaks-perf.feature
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
import {moleculeTable} from '../../bindings/footguns.js';
import {bigMolecules, chainTable, columnCreatedWithin, eachMapping, editorRounds, editorWarmed, finishedWithin, firstDrawAtMost, gridDraws, madeCurrentOnce, madeCurrentRounds, mainMemoryGrew, mapMore, mappedSome, memoryGrewAtMost, noFrameOver, noMoreDetached, noTaskOverMs, notePanel, noteSketcher, openCorpus, paneAsNoted, runOnMolecules, schemesAsNoted, schemesBounded, scrollDown, scrollTable, sizeColumn, sketcherObjectsAsNoted, twiceMolecules, watchesFromStart} from '../../bindings/leaks.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Leaks and speed: what repeated use and pathological input must never do (Crux mode)", () => {
  const session = feature(test, "features/leaks/leaks-perf.feature", import.meta.url);
  test("Opening and closing the molecule cell editor 50 times leaves Crux where it started", {tag: ["@footguns", "@batch-2", "@serial", "@catalogue-PERF-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the page records the WebAssembly memories and Crux sketchers it makes", () => watchesFromStart(page));
    await session.step(26, "And user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(32, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"c1ccc2c(c1)cc(C(=O)NCC)n2\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "c1ccc2c(c1)cc(C(=O)NCC)n2"));
    await session.step(33, "And user has opened and closed the cell editor of \"cell 1 of molecule\" 5 times", () => editorWarmed(page, "cell 1 of molecule", 5));
    await session.step(34, "And user notes the live Crux sketcher objects and the wasm memory of every Crux instance", () => noteSketcher(page));
    await session.step(35, "When user opens and closes the cell editor of \"cell 1 of molecule\" 50 times, drawing a bond each time", () => editorRounds(page, "cell 1 of molecule", 50));
    await session.step(36, "Then the live Crux sketcher objects are as many as noted", () => sketcherObjectsAsNoted(page));
    await session.step(37, "And the wasm memory of every Crux instance grew by at most 1 MiB", () => memoryGrewAtMost(page, 1));
    await session.step(38, "And the page holds no more detached sketcher roots than before", () => noMoreDetached(page));
  });
  test("Switching the current reaction cell 50 times leaves the Reaction panel's objects where they started", {tag: ["@footguns", "@batch-2", "@serial", "@catalogue-PERF-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the page records the WebAssembly memories and Crux sketchers it makes", () => watchesFromStart(page));
    await session.step(26, "And user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(42, "Given user opens the reaction corpus \"files/tests/reactions-crux.csv\" of Chem", () => openCorpus(page, "files/tests/reactions-crux.csv"));
    await session.step(43, "And user has made rows 1 to 10 current in turn", () => madeCurrentOnce(page, 1, 10));
    await session.step(44, "And user notes the live reaction schemes, the Reaction panel's canvases and their size listeners", () => notePanel(page));
    await session.step(45, "When user makes rows 1 to 10 current in turn, 5 times", () => madeCurrentRounds(page, 1, 10, 5));
    await session.step(47, "Then the live reaction schemes are as many as noted", () => schemesAsNoted(page));
    await session.step(48, "And the Reaction panel's canvases and size listeners are as many as noted", () => paneAsNoted(page));
    await session.step(49, "And the wasm memory of the main thread's Crux grew by at most 1 MiB", () => mainMemoryGrew(page, 1));
  });
  test("A 2,000-row reaction column scrolls without a frame over 100 ms", {tag: ["@footguns", "@batch-2", "@serial", "@catalogue-PERF-02"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the page records the WebAssembly memories and Crux sketchers it makes", () => watchesFromStart(page));
    await session.step(26, "And user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(53, "Given user opens a table of 2,000 reactions and routes made from the reaction corpus", () => scrollTable(page));
    await session.step(54, "And the reaction column of the grid is 600 px wide and its rows 150 px tall", () => sizeColumn(page, 600, 150));
    await session.step(55, "When user scrolls the grid from the first row to the last, one page at a time", () => scrollDown(page));
    await session.step(56, "Then no frame took more than 100 ms", () => noFrameOver(page, 100));
    await session.step(57, "And the worst frame is at most twice the worst frame of the same grid scrolled over a molecule column", () => twiceMolecules(page));
    await session.step(58, "And the live reaction schemes are at most the cell cache's bound", () => schemesBounded(page));
  });
  test("A 60-step route draws within half a second", {tag: ["@footguns", "@batch-2", "@serial", "@catalogue-PERF-04"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the page records the WebAssembly memories and Crux sketchers it makes", () => watchesFromStart(page));
    await session.step(26, "And user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(62, "Given user opens a table \"route\" with a reaction column \"reaction\" holding a 60-step chain \"CO>>CCO>>CCCO>>…\"", () => chainTable(page, "route", "reaction", 60, "CO>>CCO>>CCCO>>…"));
    await session.step(63, "When the grid draws \"cell 1 of reaction\"", () => gridDraws(page, "cell 1 of reaction"));
    await session.step(64, "Then the cell's first draw took at most 500 ms", () => firstDrawAtMost(page, 500));
  });
  test("Mapping 20 reactions in the Reaction panel answers each within 2 s", {tag: ["@footguns", "@batch-2", "@serial", "@catalogue-PERF-05"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the page records the WebAssembly memories and Crux sketchers it makes", () => watchesFromStart(page));
    await session.step(26, "And user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(72, "Given user opens the reaction corpus \"files/tests/reactions-crux.csv\" of Chem", () => openCorpus(page, "files/tests/reactions-crux.csv"));
    await session.step(73, "And user has mapped 5 reactions in the Reaction panel", () => mappedSome(page, 5));
    await session.step(74, "When user maps 20 more reactions in the Reaction panel, one after another", () => mapMore(page, 20));
    await session.step(75, "Then each mapping answered within 2 s", () => eachMapping(page, 2));
  });
  test("Transform on 50,000 molecules keeps the tab responsive and reports its counts", {tag: ["@footguns", "@batch-2", "@serial", "@catalogue-PERF-06"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given the page records the WebAssembly memories and Crux sketchers it makes", () => watchesFromStart(page));
    await session.step(26, "And user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(79, "Given user opens a table of 50,000 molecules \"smiles_200K\" cut to 50,000 rows", () => bigMolecules(page, "smiles_200K"));
    await session.step(80, "When user runs the transformation \"Aromatic Bromination\" on the molecule column", () => runOnMolecules(page, "Aromatic Bromination"));
    await session.step(83, "Then no main-thread task took more than 200 ms while it ran, the products column's creation aside", () => noTaskOverMs(page, 200));
    await session.step(84, "And the platform created the products column in one task of at most 500 ms", () => columnCreatedWithin(page, 500));
    await session.step(85, "And it finished within 60 s with a message counting products, errors and rows with no reaction", () => finishedWithin(page, 60));
  });
});
