/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/markup/markup-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.markup]
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
import {hasNoColumn, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {deleteSelected, removeColumn, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingIs, readingReads, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Markup edge cases — HTML in a cell value, a referenced column removed, every row deleted", () => {
  const session = feature(test, "features/viewers/markup/markup-edge-cases.feature", import.meta.url);
  test("HTML in a substituted cell value is shown as text and renders nothing", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "Given user opens a table \"markup-html\" with:", () => openTableOf(page, "markup-html", [["NAME","SCORE"],["<b>bold</b> <i>it</i>","1"]]), [["NAME","SCORE"],["<b>bold</b> <i>it</i>","1"]]);
    await session.step(19, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Name: ${NAME}"]]), [["content","Name: ${NAME}"]]);
    await session.step(21, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(22, "Then the \"current row\" reading of markup viewer should be 1", () => readingIs(page, "current row", el("markup viewer"), 1));
    await session.step(23, "And the \"text\" reading of markup viewer should be \"Name: <b>bold</b> <i>it</i>\"", () => readingReads(page, "text", el("markup viewer"), "Name: <b>bold</b> <i>it</i>"));
    await session.step(24, "And the \"bold\" reading of markup viewer should be 0", () => readingIs(page, "bold", el("markup viewer"), 0));
    await session.step(25, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column removed while the content references it is left as written", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: ${AGE} Sex: ${SEX}"]]), [["content","Age: ${AGE} Sex: ${SEX}"]]);
    await session.step(32, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(33, "Then the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
    await session.step(34, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(35, "Then the table should not have a column \"SEX\"", () => hasNoColumn(page, "SEX"));
    await session.step(36, "And the \"text\" reading of markup viewer should be \"Age: 26 Sex: ${SEX}\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: ${SEX}"));
    await session.step(37, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Once every row of the table is deleted the viewer has no current row and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(42, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Rows: #{t.rowCount}"]]), [["content","Rows: #{t.rowCount}"]]);
    await session.step(44, "Then the \"text\" reading of markup viewer should be \"Rows: 1000\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000"));
    await session.step(45, "When user selects all rows", () => selectAllRows(page));
    await session.step(46, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(47, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(48, "And the \"text\" reading of markup viewer should be \"Rows: 0\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 0"));
    await session.step(49, "And the \"current row\" reading of markup viewer should be 0", () => readingIs(page, "current row", el("markup viewer"), 0));
    await session.step(50, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
});
