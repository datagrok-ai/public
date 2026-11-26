/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/07-r-group-analysis.feature
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
import {clickOn, finishedUpdating, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnMatching, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/07-r-group-analysis.feature", import.meta.url);
  test("Analyze R-groups", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens sar_small dataset", () => openDataset(page, ds("sar_small")));
    await session.step(20, "And the grid shows only the columns \"smiles, ID, In-vivo Activity\"", () => gridShowsColumns(page, "smiles, ID, In-vivo Activity"));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."), undefined, "Chem > Analyze > R-Groups Analysis");
    await session.step(26, "Then \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
    await session.step(28, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")), undefined, "Draw the core, or let MCS find the common one");
    await session.step(29, "Then \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
    await session.step(30, "When the video says \"MCS: the largest substructure all the molecules share\" at the top", () => videoSaysAt(page, "MCS: the largest substructure all the molecules share", "top"));
    await session.step(31, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
    await session.step(32, "Then a new column matching \"^R1\" should have been added", () => newColumnMatching(page, "^R1"));
    await session.step(33, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
    await session.step(34, "When the video says \"One column per R-group, and a plot of how they combine\"", () => videoSays(page, "One column per R-group, and a plot of how they combine"));
    await session.step(35, "Then the table should have 200 rows", () => rowCount(page, 200));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
});
