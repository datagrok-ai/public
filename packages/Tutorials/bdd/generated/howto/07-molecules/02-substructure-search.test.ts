/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/02-substructure-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, pressKeyIn, shouldBe, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/02-substructure-search.feature", import.meta.url);
  test("Find and highlight a substructure", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(24, "And the grid shows only the columns \"molregno, canonical_smiles\"", () => gridShowsColumns(page, "molregno, canonical_smiles"));
    await session.step(25, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(26, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(27, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(29, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."), undefined, "Chem > Search > Substructure Search");
    await session.step(30, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(31, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(32, "When the video says \"Draw the fragment, or type its SMILES\" at the top", () => videoSaysAt(page, "Draw the fragment, or type its SMILES", "top"));
    await session.step(33, "And user clicks on molecule input of sketcher dialog", () => clickOn(page, el("molecule input of sketcher dialog")));
    await session.step(34, "And user types \"C(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "C(=O)O", el("molecule input of sketcher dialog")));
    await session.step(35, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(36, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(37, "Then 314 rows should pass the filter", () => filterPasses(page, 314));
    await session.step(38, "And the \"structure of canonical_smiles\" reading of filter panel should be \"C(=O)O\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "C(=O)O"));
    await session.step(39, "When the video says \"314 of 1000 molecules contain it; the match is drawn in each\"", () => videoSays(page, "314 of 1000 molecules contain it; the match is drawn in each"));
    await session.step(41, "And user hovers over the \"card canonical_smiles\" area of filter panel", () => hoverArea(page, "card canonical_smiles", el("filter panel")), undefined, "The fragment stays on the filter panel; click it to change it");
    await session.step(43, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")), undefined, "Reset brings every molecule back");
    await session.step(44, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(45, "When the video says \"To mark a fragment without hiding any rows, highlight it\"", () => videoSays(page, "To mark a fragment without hiding any rows, highlight it"));
    await session.step(47, "And user clicks on the \"header canonical_smiles\" area of grid", () => clickArea(page, "header canonical_smiles", el("grid")), undefined, "Click the column's header, then open the context panel");
    await session.step(48, "And user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")));
    await session.step(49, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
    await session.step(51, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")), undefined, "Chemistry > Highlight > Sketch");
    await session.step(52, "And user expands Highlight accordion header in context panel", () => expand(page, el("Highlight accordion header in context panel")));
    await session.step(53, "And user clicks on \"Sketch\" text in \"Highlight\" pane in context panel", () => clickOn(page, el("\"Sketch\" text in \"Highlight\" pane in context panel")));
    await session.step(54, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(56, "When user clicks on molecule input of sketcher dialog", () => clickOn(page, el("molecule input of sketcher dialog")), undefined, "Here a pyridine ring");
    await session.step(57, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
    await session.step(58, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(59, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(60, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(61, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(62, "When the video says \"The ring is colored wherever it occurs; all 1000 rows stay\"", () => videoSays(page, "The ring is colored wherever it occurs; all 1000 rows stay"));
    await session.step(63, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
});
