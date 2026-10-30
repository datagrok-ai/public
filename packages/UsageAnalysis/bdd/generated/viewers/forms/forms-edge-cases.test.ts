/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms, public#4452]
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
import {clearSelection, deleteSelected, filterPassesFewer, filterTo, removeColumn, resetFilter, rowCount, selectAllRows, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, loadLayout, noErrors, pickFromAreaContextMenu, readingAtLeast, readingIs, readingReads, saveLayout, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer edge cases", () => {
  const session = feature(test, "features/viewers/forms/forms-edge-cases.feature", import.meta.url);
  test("A table with every row deleted leaves no record card and logs no error", {tag: ["@viewers", "@realizes:viewers.forms", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(14, "Then the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await session.step(18, "Given user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(19, "Then the \"records shown\" reading of forms viewer should be at least 5", () => readingAtLeast(page, "records shown", el("forms viewer"), 5));
    await session.step(20, "When user selects all rows", () => selectAllRows(page));
    await session.step(21, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(22, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(23, "And the \"records shown\" reading of forms viewer should be 0", () => readingIs(page, "records shown", el("forms viewer"), 0));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
  test("A pinned row that is deleted leaves the pinned pane", {tag: ["@viewers", "@realizes:viewers.forms", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(14, "Then the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await session.step(28, "Given user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(29, "And user sets \"showMouseOverRow\" property of forms viewer to \"false\"", () => setProperty(page, "showMouseOverRow", el("forms viewer"), "false"));
    await session.step(30, "When user picks \"Pin Row\" from the context menu of the \"field USUBJID of card 3\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field USUBJID of card 3", el("forms viewer")));
    await session.step(31, "Then the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
    await session.step(32, "And the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
    await session.step(33, "When user clears the row selection", () => clearSelection(page));
    await session.step(34, "And user selects rows where \"USUBJID\" is \"X0273T29012500105\"", () => selectWhereIs(page, "USUBJID", "X0273T29012500105"));
    await session.step(35, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(36, "Then the table should have 999 rows", () => rowCount(page, 999));
    await session.step(37, "And the \"pinned records\" reading of forms viewer should be 0", () => readingIs(page, "pinned records", el("forms viewer"), 0));
    await session.step(38, "And the \"pinned pane shown\" reading of forms viewer should be \"false\"", () => readingReads(page, "pinned pane shown", el("forms viewer"), "false"));
    await session.step(39, "And the \"pinned values\" reading of forms viewer should be \"\"", () => readingReads(page, "pinned values", el("forms viewer"), ""));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("A pinned row the filter hides stays pinned", {tag: ["@viewers", "@realizes:viewers.forms", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(14, "Then the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await session.step(44, "Given user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(45, "And user sets \"showMouseOverRow\" property of forms viewer to \"false\"", () => setProperty(page, "showMouseOverRow", el("forms viewer"), "false"));
    await session.step(46, "When user picks \"Pin Row\" from the context menu of the \"field USUBJID of card 3\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field USUBJID of card 3", el("forms viewer")));
    await session.step(47, "Then the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
    await session.step(48, "When user clears the row selection", () => clearSelection(page));
    await session.step(49, "And user filters rows where \"SEVERITY\" is \"Low\"", () => filterTo(page, "SEVERITY", "Low"));
    await session.step(50, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(51, "And the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
    await session.step(52, "And the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
    await session.step(53, "When user resets the filter", () => resetFilter(page));
    await session.step(54, "Then the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a field column the table no longer has shows the other fields", {tag: ["@viewers", "@realizes:viewers.forms", "@edge-cases", "@known-failure:66", "@realizes:public#4452"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(14, "Then the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await knownFailure(async () => {
      await session.step(61, "When user sets \"fields\" property of forms viewer to \"AGE, HEIGHT, SEX\"", () => setProperty(page, "fields", el("forms viewer"), "AGE, HEIGHT, SEX"));
      await session.step(62, "Then the \"fields shown\" reading of forms viewer should be 3", () => readingIs(page, "fields shown", el("forms viewer"), 3));
      await session.step(63, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(64, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
      await session.step(65, "And user loads the saved layout", () => loadLayout(page));
      await session.step(66, "Then the \"fields shown\" reading of forms viewer should be 2", () => readingIs(page, "fields shown", el("forms viewer"), 2));
      await session.step(67, "And the \"fields\" reading of forms viewer should be \"AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "AGE, SEX"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    }, 66);
  });
});
