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
import {clickOn, enterInto, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn, makeLastRowCurrent, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {deleteSelected, filterPasses, filterTo, removeColumn, renameColumn, rowCount, selectAllRows, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, pickFromContextMenu, propertyShouldBe, readingIs, readingReads, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

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
  test("Edit content with new text typed and then CANCEL leaves the content as it was [dismissal=CANCEL, dismiss=user clicks on \"CANCEL\" button in Edit dialog]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(54, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: ${AGE}"]]), [["content","Age: ${AGE}"]]);
    await session.step(57, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(58, "Then the \"text\" reading of markup viewer should be \"Age: 26\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26"));
    await session.step(59, "When user picks \"Edit content...\" from the context menu of markup viewer", () => pickFromContextMenu(page, "Edit content...", el("markup viewer")));
    await session.step(60, "Then Edit dialog should be visible", () => shouldBe(page, el("Edit dialog"), "visible"));
    await session.step(61, "When user enters \"Weight: ${WEIGHT}\" into input in Edit dialog", () => enterInto(page, "Weight: ${WEIGHT}", el("input in Edit dialog")));
    await session.step(62, "And user clicks on \"CANCEL\" button in Edit dialog", () => clickOn(page, el("\"CANCEL\" button in Edit dialog")));
    await session.step(63, "Then Edit dialog should be absent", () => shouldBe(page, el("Edit dialog"), "absent"));
    await session.step(64, "And \"content\" property of markup viewer should be \"Age: ${AGE}\"", () => propertyShouldBe(page, "content", el("markup viewer"), "Age: ${AGE}"));
    await session.step(65, "And the \"text\" reading of markup viewer should be \"Age: 26\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("Edit content with new text typed and then Escape leaves the content as it was [dismissal=Escape, dismiss=user presses Escape]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(54, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: ${AGE}"]]), [["content","Age: ${AGE}"]]);
    await session.step(57, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(58, "Then the \"text\" reading of markup viewer should be \"Age: 26\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26"));
    await session.step(59, "When user picks \"Edit content...\" from the context menu of markup viewer", () => pickFromContextMenu(page, "Edit content...", el("markup viewer")));
    await session.step(60, "Then Edit dialog should be visible", () => shouldBe(page, el("Edit dialog"), "visible"));
    await session.step(61, "When user enters \"Weight: ${WEIGHT}\" into input in Edit dialog", () => enterInto(page, "Weight: ${WEIGHT}", el("input in Edit dialog")));
    await session.step(62, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(63, "Then Edit dialog should be absent", () => shouldBe(page, el("Edit dialog"), "absent"));
    await session.step(64, "And \"content\" property of markup viewer should be \"Age: ${AGE}\"", () => propertyShouldBe(page, "content", el("markup viewer"), "Age: ${AGE}"));
    await session.step(65, "And the \"text\" reading of markup viewer should be \"Age: 26\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("With the current row filtered out the placeholder shows empty", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup", "@known-failure:82", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(75, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(76, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: [${AGE}]"]]), [["content","Age: [${AGE}]"]]);
      await session.step(78, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(79, "Then the \"text\" reading of markup viewer should be \"Age: [26]\"", () => readingReads(page, "text", el("markup viewer"), "Age: [26]"));
      await session.step(80, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(81, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(82, "And the \"text\" reading of markup viewer should be \"Age: []\"", () => readingReads(page, "text", el("markup viewer"), "Age: []"));
      await session.step(83, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    }, 82);
  });
  test("Once the current row is deleted the table has no current row and the placeholder shows empty", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup", "@known-failure:100", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(88, "Given user opens a table \"markup-no-current\" with:", () => openTableOf(page, "markup-no-current", [["NAME","SCORE"],["first","1"],["last","2"]]), [["NAME","SCORE"],["first","1"],["last","2"]]);
      await session.step(92, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Name: [${NAME}]"]]), [["content","Name: [${NAME}]"]]);
      await session.step(94, "And user makes the last row current", () => makeLastRowCurrent(page));
      await session.step(95, "Then the \"text\" reading of markup viewer should be \"Name: [last]\"", () => readingReads(page, "text", el("markup viewer"), "Name: [last]"));
      await session.step(96, "When user selects rows where \"NAME\" is \"last\"", () => selectWhereIs(page, "NAME", "last"));
      await session.step(97, "And user deletes the selected rows", () => deleteSelected(page));
      await session.step(98, "Then the table should have 1 row", () => rowCount(page, 1));
      await session.step(99, "And the \"current row\" reading of markup viewer should be 0", () => readingIs(page, "current row", el("markup viewer"), 0));
      await session.step(100, "And the \"text\" reading of markup viewer should be \"Name: []\"", () => readingReads(page, "text", el("markup viewer"), "Name: []"));
      await session.step(101, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    }, 100);
  });
  test("A referenced column renamed is followed by the content", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.markup", "@known-failure:113", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(106, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(107, "When user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: ${AGE}"]]), [["content","Age: ${AGE}"]]);
      await session.step(109, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(110, "Then the \"text\" reading of markup viewer should be \"Age: 26\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26"));
      await session.step(111, "When user renames \"AGE\" column to \"YEARS\"", () => renameColumn(page, "AGE", "YEARS"));
      await session.step(112, "And user makes row 4 current", () => makeRowCurrent(page, 4));
      await session.step(113, "Then the \"text\" reading of markup viewer should be \"Age: 45\"", () => readingReads(page, "text", el("markup viewer"), "Age: 45"));
      await session.step(114, "And \"content\" property of markup viewer should be \"Age: ${YEARS}\"", () => propertyShouldBe(page, "content", el("markup viewer"), "Age: ${YEARS}"));
      await session.step(115, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    }, 113);
  });
});
