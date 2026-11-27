/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/07-substructure-filter.feature
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
import {clickOn, pressKeyIn, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openCardSettings, pickSearchType} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/07-substructure-filter.feature", import.meta.url);
  test("Filter molecules by a substructure", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(16, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(23, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(25, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")), undefined, "Sketch opens the molecule editor");
    await session.step(26, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(28, "When user clicks on molecule input of sketcher dialog", () => clickOn(page, el("molecule input of sketcher dialog")), undefined, "Draw a structure, or type or paste its SMILES");
    await session.step(29, "And user types \"C1CCNC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCNC1", el("molecule input of sketcher dialog")));
    await session.step(30, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(31, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(32, "Then 33 rows should pass the filter", () => filterPasses(page, 33));
    await session.step(33, "And the \"structure of Structure\" reading of filter panel should be \"C1CCNC1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "C1CCNC1"));
    await session.step(34, "When the video says \"33 of 100 molecules contain the ring; it is highlighted in each\"", () => videoSays(page, "33 of 100 molecules contain the ring; it is highlighted in each"));
    await session.step(36, "And user opens the settings of the \"Structure\" filter card", () => openCardSettings(page, "Structure"), undefined, "The search type can flip it to Not contains");
    await session.step(37, "And user picks search type \"Not contains\" in the \"Structure\" filter card", () => pickSearchType(page, "Not contains", "Structure"));
    await session.step(38, "Then the \"search type of Structure\" reading of filter panel should be \"Not contains\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Not contains"));
    await session.step(39, "And 67 rows should pass the filter", () => filterPasses(page, 67));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
});
