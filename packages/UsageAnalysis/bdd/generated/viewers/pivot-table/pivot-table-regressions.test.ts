/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table, public#2606, public#3470]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCategoricalFilter, selectWhereIs, setLinkKeyPair, setLinkTables, tableFilterCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noBalloons, noErrors, pickFromAreaContextMenu, readingDoesNotRead, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table regressions from closed issues", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-regressions.feature", import.meta.url);
  test("At Row Source Filtered the inner grid keeps its colouring and its hidden column through a filter (public#2606)", {tag: ["@viewers", "@realizes:viewers.pivot-table", "@realizes:public#2606"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(14, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(19, "When user sets \"Row Source\" property of pivot table viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Filtered"));
    await session.step(20, "Then the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
    await session.step(21, "When user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None avg(AGE)", el("pivot table viewer")));
    await session.step(22, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
    await session.step(23, "When user picks \"Grid > Hide\" from the context menu of the \"grid header High avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Hide", "grid header High avg(AGE)", el("pivot table viewer")));
    await session.step(24, "Then the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
    await session.step(25, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
    await session.step(26, "Then the \"rows shown\" reading of pivot table viewer should be 447", () => readingIs(page, "rows shown", el("pivot table viewer"), 447));
    await session.step(27, "And the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
    await session.step(28, "And the \"color of grid cell 2 of None avg(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
    await session.step(29, "And the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("At Row Source Selected the inner grid keeps its colouring and its hidden column through a selection (public#2606)", {tag: ["@viewers", "@realizes:viewers.pivot-table", "@realizes:public#2606"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(14, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(35, "Given user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
    await session.step(36, "When user sets \"Row Source\" property of pivot table viewer to \"Selected\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Selected"));
    await session.step(37, "Then the \"rows shown\" reading of pivot table viewer should be 553", () => readingIs(page, "rows shown", el("pivot table viewer"), 553));
    await session.step(38, "When user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None avg(AGE)", el("pivot table viewer")));
    await session.step(39, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
    await session.step(40, "When user picks \"Grid > Hide\" from the context menu of the \"grid header High avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Hide", "grid header High avg(AGE)", el("pivot table viewer")));
    await session.step(41, "Then the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
    await session.step(42, "When user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
    await session.step(43, "Then the \"rows shown\" reading of pivot table viewer should be 447", () => readingIs(page, "rows shown", el("pivot table viewer"), 447));
    await session.step(44, "And the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
    await session.step(45, "And the \"color of grid cell 2 of None avg(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
    await session.step(46, "And the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("Link Tables links the source with the table the pivot published (public#3470)", {tag: ["@viewers", "@realizes:viewers.pivot-table", "@realizes:public#3470"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(14, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(52, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
    await session.step(53, "Then table \"demog-1000 aggregation\" should be open", () => tableOpen(page, "demog-1000 aggregation"));
    await session.step(54, "When user picks \"Data > Link Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Link Tables..."));
    await session.step(55, "And user sets the tables of the Link Tables dialog to \"demog-1000\" and \"demog-1000 aggregation\"", () => setLinkTables(page, "demog-1000", "demog-1000 aggregation"));
    await session.step(56, "And user sets key columns 1 of the Link Tables dialog to \"DIS_POP\" and \"DIS_POP\"", () => setLinkKeyPair(page, 1, "DIS_POP", "DIS_POP"));
    await session.step(57, "And user clicks on LINK button in \"Link Tables\" dialog", () => clickOn(page, el("LINK button in \"Link Tables\" dialog")));
    await session.step(58, "Then UNLINK button in \"Link Tables\" dialog should be visible", () => shouldBe(page, el("UNLINK button in \"Link Tables\" dialog"), "visible"));
    await session.step(59, "And 1 rows of table \"demog-1000 aggregation\" should pass the filter", () => tableFilterCount(page, 1, "demog-1000 aggregation"));
    await session.step(60, "When user clicks on CLOSE button in \"Link Tables\" dialog", () => clickOn(page, el("CLOSE button in \"Link Tables\" dialog")));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
