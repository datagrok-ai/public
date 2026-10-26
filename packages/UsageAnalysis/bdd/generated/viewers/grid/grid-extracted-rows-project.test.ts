/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-extracted-rows-project.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
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
import {pressKey, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectFirstRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, dialogCloses, openDataset, openProjectWithTable, openSaveDialog, ownProjectGone, projectsOnServer, saveFromDialog, simpleModeOff, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft, openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, readingReads, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A project holding a table of extracted rows reopens", () => {
  const session = feature(test, "features/viewers/grid/grid-extracted-rows-project.feature", import.meta.url);
  test("Extracted rows saved in a project come back with their grid and no error", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the user's own project \"BDD-Extracted-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Extracted-{time}")));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "When user selects the first 5 rows", () => selectFirstRows(page, 5));
    await session.step(20, "And user picks \"Select > Extract Selected Rows\" from the top menu", () => pickFromTopMenu(page, "Select > Extract Selected Rows"));
    await session.step(21, "Then the open table views should be exactly \"demog-1000, demog-1000 (2)\"", () => openTableViewsExactly(page, "demog-1000, demog-1000 (2)"));
    await session.step(22, "And grid should show 5 rows", () => showsRows(page, el("grid"), 5));
    await session.step(23, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
    await session.step(24, "And user types \"BDD-Extracted-{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDD-Extracted-{time}"), el("text input in \"Save project\" dialog")));
    await session.step(25, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(26, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(27, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(28, "Then 1 project named \"BDD-Extracted-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Extracted-{time}")));
    await session.step(29, "When user closes all views", () => closeAllViews(page));
    await session.step(30, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(31, "When user opens the \"BDD-Extracted-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Extracted-{time}")));
    await session.step(32, "Then the open table views should be exactly \"demog-1000, demog-1000 (2)\"", () => openTableViewsExactly(page, "demog-1000, demog-1000 (2)"));
    await session.step(33, "When user switches to the \"demog-1000 (2)\" table view", () => switchTableView(page, "demog-1000 (2)"));
    await session.step(34, "Then grid should show 5 rows", () => showsRows(page, el("grid"), 5));
    await session.step(35, "And the \"text of cell 1 of USUBJID\" reading of grid should be \"X0273T21000300003\"", () => readingReads(page, "text of cell 1 of USUBJID", el("grid"), "X0273T21000300003"));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
