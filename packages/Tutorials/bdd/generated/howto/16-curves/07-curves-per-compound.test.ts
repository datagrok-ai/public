/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/07-curves-per-compound.feature
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
import {clickOn, expand, hoverOver, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {allOfFiltered, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, closeContextMenu, hasArea, noErrors, pickFromAreaContextMenu, readingAtLeast, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInAreaSelector, pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/07-curves-per-compound.feature", import.meta.url);
  test("Curves per compound in a trellis plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(23, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(24, "When the video says \"A trellis plot draws one small chart per group\"", () => videoSays(page, "A trellis plot draws one small chart per group"));
    await session.step(25, "And user clicks on trellis plot icon on toolbox", () => clickOn(page, el("trellis plot icon on toolbox")));
    await session.step(26, "Then trellis plot viewer should be added to the open tableview", () => viewerAdded(page, "trellis plot"));
    await session.step(28, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(29, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(31, "When user picks \"Curves\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Curves", el("trellis plot viewer")), undefined, "Curves: the chart drawn in every cell");
    await session.step(32, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"MultiCurveViewer\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "MultiCurveViewer"));
    await session.step(34, "When user picks \"Compound ID\" in the column selector at the \"x selector 1\" area of trellis plot viewer", () => pickInAreaSelector(page, "Compound ID", "x selector 1", el("trellis plot viewer")), undefined, "Split by compound across");
    await session.step(36, "And user picks \"Reset Y columns\" from the context menu of the \"y selector 1\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Reset Y columns", "y selector 1", el("trellis plot viewer")), undefined, "Drop the split down: right-click, Reset Y columns");
    await session.step(37, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(38, "Then the \"cells\" reading of trellis plot viewer should be at least 20", () => readingAtLeast(page, "cells", el("trellis plot viewer"), 20));
    await session.step(39, "And trellis plot viewer should have a \"cell compound_id 20\" area", () => hasArea(page, el("trellis plot viewer"), "cell compound_id 20"));
    await session.step(40, "When the video says \"Each compound's runs, overlaid in a cell of its own\"", () => videoSays(page, "Each compound's runs, overlaid in a cell of its own"));
    await session.step(42, "When user hovers over trellis plot viewer", () => hoverOver(page, el("trellis plot viewer")), undefined, "Settings: the gear on the plot's title bar");
    await session.step(43, "And user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
    await session.step(45, "And user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc > On Click: Filter");
    await session.step(46, "And user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
    await session.step(48, "And user clicks on the \"cell body compound_id 20\" area of trellis plot viewer", () => clickArea(page, "cell body compound_id 20", el("trellis plot viewer")), undefined, "Click a compound's chart: the table shows only its runs");
    await session.step(49, "Then 4 rows should pass the filter", () => filterPasses(page, 4));
    await session.step(50, "And all rows where \"Compound ID\" is \"compound_id 20\" should pass the filter", () => allOfFiltered(page, "Compound ID", "compound_id 20"));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
});
