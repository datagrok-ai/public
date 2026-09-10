/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/markup/markup-row-binding.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.markup]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {currentRowValue, filterPasses, selectAllRows, selectNoRows, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingIs, readingReads, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingExcludes, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Markup column references, table expressions and the current row", () => {
  const session = feature(test, "features/viewers/markup/markup-row-binding.feature", import.meta.url);
  test("Markup column references, table expressions and the current row", {tag: ["@journey", "@viewers", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: ${AGE} Sex: ${SEX}"]]));
    await session.step(23, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(24, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(25, "And the \"current row\" reading of markup viewer should be 1", () => readingIs(page, "current row", el("markup viewer"), 1));
    await session.step(26, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await run.scenario("The references render the current row's values and follow it", async () => {
      await session.step(29, "Then \"AGE\" of the current row should be \"26\"", () => currentRowValue(page, "AGE", "26"));
      await session.step(30, "And the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
      await session.step(31, "When user makes row 4 current", () => makeRowCurrent(page, 4));
      await session.step(32, "Then the \"current row\" reading of markup viewer should be 4", () => readingIs(page, "current row", el("markup viewer"), 4));
      await session.step(33, "And \"AGE\" of the current row should be \"45\"", () => currentRowValue(page, "AGE", "45"));
      await session.step(34, "And the \"text\" reading of markup viewer should be \"Age: 45 Sex: M\"", () => readingReads(page, "text", el("markup viewer"), "Age: 45 Sex: M"));
      await session.step(35, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(36, "Then the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A reference to a column that does not exist is left exactly as written", async () => {
      await session.step(40, "When user sets \"content\" property of markup viewer to \"Known: ${AGE} Unknown: ${NO_SUCH_COLUMN}\"", () => setProperty(page, "content", el("markup viewer"), "Known: ${AGE} Unknown: ${NO_SUCH_COLUMN}"));
      await session.step(41, "Then the \"text\" reading of markup viewer should be \"Known: 26 Unknown: ${NO_SUCH_COLUMN}\"", () => readingReads(page, "text", el("markup viewer"), "Known: 26 Unknown: ${NO_SUCH_COLUMN}"));
      await session.step(42, "And the \"text\" reading of markup viewer should not include the text \"${AGE}\"", () => readingExcludes(page, "text", el("markup viewer"), "${AGE}"));
      await session.step(43, "When user sets \"content\" property of markup viewer to \"Age: ${AGE} Sex: ${SEX}\"", () => setProperty(page, "content", el("markup viewer"), "Age: ${AGE} Sex: ${SEX}"));
      await session.step(44, "Then the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The table expressions count the rows and the selection, and follow both", async () => {
      await session.step(48, "When user sets \"content\" property of markup viewer to \"Rows: #{t.rowCount} Selected: #{t.selection.trueCount}\"", () => setProperty(page, "content", el("markup viewer"), "Rows: #{t.rowCount} Selected: #{t.selection.trueCount}"));
      await session.step(49, "Then the \"text\" reading of markup viewer should be \"Rows: 1000 Selected: 0\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000 Selected: 0"));
      await session.step(50, "When user selects all rows", () => selectAllRows(page));
      await session.step(51, "Then 1000 rows should be selected", () => selectedRowCount(page, 1000));
      await session.step(52, "And the \"text\" reading of markup viewer should be \"Rows: 1000 Selected: 1000\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000 Selected: 1000"));
      await session.step(53, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(54, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(55, "And the \"text\" reading of markup viewer should be \"Rows: 1000 Selected: 553\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000 Selected: 553"));
      await session.step(56, "When user selects no rows", () => selectNoRows(page));
      await session.step(57, "Then the \"text\" reading of markup viewer should be \"Rows: 1000 Selected: 0\"", () => readingReads(page, "text", el("markup viewer"), "Rows: 1000 Selected: 0"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A heading and a list are rendered around the substituted values", async () => {
      await session.step(61, "When user sets \"content\" property of markup viewer to \"# Demographics\\n\\n* Age: ${AGE}\\n* Sex: ${SEX}\"", () => setProperty(page, "content", el("markup viewer"), "# Demographics\\n\\n* Age: ${AGE}\\n* Sex: ${SEX}"));
      await session.step(62, "Then the \"mode\" reading of markup viewer should be \"Markup\"", () => readingReads(page, "mode", el("markup viewer"), "Markup"));
      await session.step(63, "And the \"heading 1\" reading of markup viewer should be \"Demographics\"", () => readingReads(page, "heading 1", el("markup viewer"), "Demographics"));
      await session.step(64, "And the \"list items\" reading of markup viewer should be 2", () => readingIs(page, "list items", el("markup viewer"), 2));
      await session.step(65, "And the \"text\" reading of markup viewer should include the text \"Age: 26\"", () => readingIncludes(page, "text", el("markup viewer"), "Age: 26"));
      await session.step(66, "And the \"text\" reading of markup viewer should include the text \"Sex: F\"", () => readingIncludes(page, "text", el("markup viewer"), "Sex: F"));
      await session.step(67, "When user sets \"content\" property of markup viewer to \"Age: ${AGE} Sex: ${SEX}\"", () => setProperty(page, "content", el("markup viewer"), "Age: ${AGE} Sex: ${SEX}"));
      await session.step(68, "Then the \"heading 1\" reading of markup viewer should be \"\"", () => readingReads(page, "heading 1", el("markup viewer"), ""));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
