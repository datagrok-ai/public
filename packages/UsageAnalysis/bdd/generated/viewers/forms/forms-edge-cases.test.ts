/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
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
import {clearSelection, deleteSelected, rowCount, selectAllRows, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, pickFromAreaContextMenu, readingAtLeast, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer edge cases", () => {
  const session = feature(test, "features/viewers/forms/forms-edge-cases.feature", import.meta.url);
  test("A table with every row deleted leaves no record card and logs no error", {tag: ["@viewers", "@realizes:viewers.forms", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(13, "Then the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await session.step(17, "Given user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(18, "Then the \"records shown\" reading of forms viewer should be at least 5", () => readingAtLeast(page, "records shown", el("forms viewer"), 5));
    await session.step(19, "When user selects all rows", () => selectAllRows(page));
    await session.step(20, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(21, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(22, "And the \"records shown\" reading of forms viewer should be 0", () => readingIs(page, "records shown", el("forms viewer"), 0));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
  test("A pinned row that is deleted leaves the pinned pane", {tag: ["@viewers", "@realizes:viewers.forms", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(13, "Then the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await session.step(27, "Given user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(28, "And user sets \"showMouseOverRow\" property of forms viewer to \"false\"", () => setProperty(page, "showMouseOverRow", el("forms viewer"), "false"));
    await session.step(29, "When user picks \"Pin Row\" from the context menu of the \"field USUBJID of card 3\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field USUBJID of card 3", el("forms viewer")));
    await session.step(30, "Then the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
    await session.step(31, "And the \"pinned values\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "pinned values", el("forms viewer"), "X0273T29012500105"));
    await session.step(32, "When user clears the row selection", () => clearSelection(page));
    await session.step(33, "And user selects rows where \"USUBJID\" is \"X0273T29012500105\"", () => selectWhereIs(page, "USUBJID", "X0273T29012500105"));
    await session.step(34, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(35, "Then the table should have 999 rows", () => rowCount(page, 999));
    await session.step(36, "And the \"pinned records\" reading of forms viewer should be 0", () => readingIs(page, "pinned records", el("forms viewer"), 0));
    await session.step(37, "And the \"pinned pane shown\" reading of forms viewer should be \"false\"", () => readingReads(page, "pinned pane shown", el("forms viewer"), "false"));
    await session.step(38, "And the \"pinned values\" reading of forms viewer should be \"\"", () => readingReads(page, "pinned values", el("forms viewer"), ""));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
});
