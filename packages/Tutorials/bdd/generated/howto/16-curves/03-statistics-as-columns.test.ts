/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/16-curves/03-statistics-as-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
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
import {clickBesideCurve} from '../../../bindings/howto-curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, displayedInRow, distinctValues, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, contextPanelOpen, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, hoverArea, noErrors, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Curves", () => {
  const session = feature(test, "features/howto/16-curves/03-statistics-as-columns.feature", import.meta.url);
  test("IC50 and other statistics as columns", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the grid shows only the columns \"Molecule, Fitted Curve, Compound|Assay|Target\"", () => gridShowsColumns(page, "Molecule, Fitted Curve, Compound|Assay|Target"));
    await session.step(23, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(25, "When user clicks on the \"cell 1 of Fitted Curve\" area of grid beside its curve", () => clickBesideCurve(page, "cell 1 of Fitted Curve", el("grid")), undefined, "Click a curve, then open Fit in the context panel");
    await session.step(26, "Then \"Fit\" pane in context panel should be visible", () => shouldBe(page, el("\"Fit\" pane in context panel"), "visible"));
    await session.step(27, "When user expands Fit accordion header in context panel", () => expand(page, el("Fit accordion header in context panel")));
    await session.step(28, "Then \"IC50\" input in context panel should be visible", () => shouldBe(page, el("\"IC50\" input in context panel"), "visible"));
    await session.step(29, "When the video says \"Fit: this curve's statistics, from its fitted parameters\"", () => videoSays(page, "Fit: this curve's statistics, from its fitted parameters"));
    await session.step(31, "And user hovers over first \"Calculate ic50 for the whole column\" icon in context panel", () => hoverOver(page, el("first \"Calculate ic50 for the whole column\" icon in context panel")), undefined, "The + next to a value adds it for every curve");
    await session.step(32, "And user clicks on first \"Calculate ic50 for the whole column\" icon in context panel", () => clickOn(page, el("first \"Calculate ic50 for the whole column\" icon in context panel")));
    await session.step(33, "Then the table should have a column \"Fitted Curve 1 ic50\"", () => hasColumn(page, "Fitted Curve 1 ic50"));
    await session.step(34, "And \"Fitted Curve 1 ic50\" column should have type \"double\"", () => columnType(page, "Fitted Curve 1 ic50", "double"));
    await session.step(35, "And \"Fitted Curve 1 ic50\" column should have at least 40 distinct values", () => distinctValues(page, "Fitted Curve 1 ic50", 40));
    await session.step(36, "When the video says \"A new column with every curve's IC50, too small to read yet\"", () => videoSays(page, "A new column with every curve's IC50, too small to read yet"));
    await session.step(38, "And user clicks on first \"Calculate rSquared for the whole column\" icon in context panel", () => clickOn(page, el("first \"Calculate rSquared for the whole column\" icon in context panel")), undefined, "R² tells how well each curve fits its points");
    await session.step(39, "Then the table should have a column \"Fitted Curve 1 rSquared\"", () => hasColumn(page, "Fitted Curve 1 rSquared"));
    await session.step(41, "When user picks \"Format > scientific (1.63E-7)\" from the context menu of the \"header Fitted Curve 1 ic50\" area of grid", () => pickFromAreaContextMenu(page, "Format > scientific (1.63E-7)", "header Fitted Curve 1 ic50", el("grid")), undefined, "IC50s in molar are tiny: show them in scientific notation");
    await session.step(42, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(43, "Then the \"Fitted Curve 1 ic50\" cell of row 1 should be displayed as \"1.63E-7\"", () => displayedInRow(page, "Fitted Curve 1 ic50", 1, "1.63E-7"));
    await session.step(45, "When user picks \"Sort > Ascending\" from the context menu of the \"header Fitted Curve 1 rSquared\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Ascending", "header Fitted Curve 1 rSquared", el("grid")), undefined, "Sort by R²: the poorest fits come first");
    await session.step(46, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(47, "Then the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
    await session.step(50, "When user hovers over the \"cell 7 of Fitted Curve 1 rSquared\" area of grid", () => hoverArea(page, "cell 7 of Fitted Curve 1 rSquared", el("grid")), undefined, "R² near 0: the curve does not fit these points");
    await session.step(51, "Then no errors should have been logged", () => noErrors(page));
  });
});
