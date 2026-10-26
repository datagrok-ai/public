/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, public#1071, public#1887]
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
import {categoricalColorIs, colorCategorical, colorCodedCategorically, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openDatasetRowsAs, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {loadLayout, noErrors, pickFromAreaContextMenu, readingReads, readingsDiffer, readingsEqual, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Colour coding regressions from closed issues", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-regressions.feature", import.meta.url);
  test("After the layout comes back, Off in the header menu switches the colouring off (public#1071)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#1071"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(16, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(21, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#3366CC"],["F","#CC6699"]]), [["M","#3366CC"],["F","#CC6699"]]);
    await session.step(24, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(25, "And user closes all views", () => closeAllViews(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(28, "When user loads the saved layout", () => loadLayout(page));
    await session.step(29, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
    await session.step(30, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
    await session.step(31, "When user picks \"Color Coding > Off\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header SEX", el("grid")));
    await session.step(32, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(33, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
    await session.step(34, "And the \"color of cell 1 of SEX\" reading of grid should be \"#ffffff\"", () => readingReads(page, "color of cell 1 of SEX", el("grid"), "#ffffff"));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
  });
  test("A colouring picked up in one table applies to a column of another table (public#1887)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#1887"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(16, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(40, "Given user opens demog-1000 dataset keeping the first 100 rows as \"demog-100\"", () => openDatasetRowsAs(page, ds("demog-1000"), 100, "demog-100"));
    await session.step(41, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(42, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#3366CC"],["Other","#CC6699"]]), [["Caucasian","#3366CC"],["Other","#CC6699"]]);
    await session.step(45, "And user picks \"Pick Up / Apply > Pick Up\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Pick Up", "header RACE", el("grid")));
    await session.step(46, "And user switches to the \"demog-100\" table view", () => switchTableView(page, "demog-100"));
    await session.step(47, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(48, "When user picks \"Pick Up / Apply > Apply\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Apply", "header RACE", el("grid")));
    await session.step(49, "Then \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
    await session.step(50, "And the categorical color of \"Caucasian\" in \"RACE\" column should be \"3366CC\"", () => categoricalColorIs(page, "Caucasian", "RACE", "3366CC"));
    await session.step(51, "And the categorical color of \"Other\" in \"RACE\" column should be \"CC6699\"", () => categoricalColorIs(page, "Other", "RACE", "CC6699"));
    await session.step(52, "And the \"color of cell 1 of RACE\" reading of grid should be \"#3366cc\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#3366cc"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
});
