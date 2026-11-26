/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/04-calculate-properties.feature
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
import {check, clickOn, hoverOver, pressKeyIn, shouldBe, shouldContainText, typeInto, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueBetween} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactly, filterPasses, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/04-calculate-properties.feature", import.meta.url);
  test("Calculate properties for every molecule", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(22, "And the grid shows only the columns \"molregno, canonical_smiles\"", () => gridShowsColumns(page, "molregno, canonical_smiles"));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(27, "When user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."), undefined, "Chem > Calculate > Chemical Properties");
    await session.step(28, "Then \"Chemical Properties\" dialog should be visible", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "visible"));
    await session.step(29, "And Molecules input in \"Chemical Properties\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Chemical Properties\" dialog"), "canonical_smiles"));
    await session.step(30, "And MW input in \"Chemical Properties\" dialog should be checked", () => shouldBe(page, el("MW input in \"Chemical Properties\" dialog"), "checked"));
    await session.step(32, "When user checks \"Chemical Properties\" calculator in \"Chemical Properties\" dialog", () => check(page, el("\"Chemical Properties\" calculator in \"Chemical Properties\" dialog")), undefined, "Tick the calculator and the properties you need");
    await session.step(33, "And user checks \"Log P\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Log P\" input in \"Chemical Properties\" dialog")));
    await session.step(34, "And user checks PSA input in \"Chemical Properties\" dialog", () => check(page, el("PSA input in \"Chemical Properties\" dialog")));
    await session.step(35, "And user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
    await session.step(36, "Then \"Chemical Properties\" dialog should be hidden", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "hidden"));
    await session.step(37, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
    await session.step(38, "And a new column \"MW\" should have been added", () => newColumnNamed(page, "MW"));
    await session.step(39, "And a new column \"LogP\" should have been added", () => newColumnNamed(page, "LogP"));
    await session.step(40, "And a new column \"PSA\" should have been added", () => newColumnNamed(page, "PSA"));
    await session.step(41, "And every value of \"MW\" column should lie between 1 and 5000", () => everyValueBetween(page, "MW", 1, 5000));
    await session.step(42, "When the video says \"A new column per property, computed in your browser\" at the bottom", () => videoSaysAt(page, "A new column per property, computed in your browser", "bottom"));
    await session.step(43, "And the video says \"MW is the weight; LogP, oil versus water; PSA, the polar surface\" at the bottom", () => videoSaysAt(page, "MW is the weight; LogP, oil versus water; PSA, the polar surface", "bottom"));
    await session.step(44, "And the video says \"Most oral drugs stay under MW 500 and LogP 5\" at the bottom", () => videoSaysAt(page, "Most oral drugs stay under MW 500 and LogP 5", "bottom"));
    await session.step(46, "And user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")), undefined, "The Search box (Ctrl+F) takes a condition");
    await session.step(47, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(49, "When user types \"MW > 500\" into table search", () => typeInto(page, "MW > 500", el("table search")), undefined, "Enter keeps only the rows that match");
    await session.step(50, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(51, "Then 211 rows should pass the filter", () => filterPasses(page, 211));
    await session.step(52, "And the filter should pass exactly the rows where \"MW\" is between 500.001 and 5000", () => filterIsExactly(page, "MW", 500.001, 5000));
    await session.step(53, "When the video says \"211 of the 1000 molecules are heavier than 500\" at the bottom", () => videoSaysAt(page, "211 of the 1000 molecules are heavier than 500", "bottom"));
    await session.step(54, "And user hovers over \"Filtered: 211\" text in status bar", () => hoverOver(page, el("\"Filtered: 211\" text in status bar")));
    await session.step(56, "And user double-clicks on the \"header MW\" area of grid", () => doubleClickArea(page, "header MW", el("grid")), undefined, "Double-click a header to sort by it: heaviest first");
    await session.step(57, "Then the \"sort column\" reading of grid should be \"MW\"", () => readingReads(page, "sort column", el("grid"), "MW"));
    await session.step(58, "And the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
    await session.step(59, "When the video says \"The heaviest molecules now come first\" at the bottom", () => videoSaysAt(page, "The heaviest molecules now come first", "bottom"));
    await session.step(60, "Then 211 rows should pass the filter", () => filterPasses(page, 211));
    await session.step(61, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
});
