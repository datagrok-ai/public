/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/markup/markup-person-paths.feature
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
import {clickOn, enterInto, pressKey, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors, pickFromContextMenu, propertyShouldBe, readingIs, readingReads, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Markup — the keyboard in the grid moves it, Show Title hides its title, Edit content applies", () => {
  const session = feature(test, "features/viewers/markup/markup-person-paths.feature", import.meta.url);
  test("Markup — the keyboard in the grid moves it, Show Title hides its title, Edit content applies", {tag: ["@journey", "@viewers", "@realizes:viewers.markup"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","Age: ${AGE} Sex: ${SEX}"]]), [["content","Age: ${AGE} Sex: ${SEX}"]]);
    await session.step(19, "Then markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await run.scenario("ArrowDown in the grid moves the current row, and the markup follows it", async () => {
      await session.step(22, "When user clicks on the \"cell 1 of AGE\" area of grid", () => clickArea(page, "cell 1 of AGE", el("grid")));
      await session.step(23, "Then the \"current row\" reading of markup viewer should be 1", () => readingIs(page, "current row", el("markup viewer"), 1));
      await session.step(24, "And the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
      await session.step(25, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
      await session.step(26, "Then the \"current row\" reading of markup viewer should be 2", () => readingIs(page, "current row", el("markup viewer"), 2));
      await session.step(27, "And the \"text\" reading of markup viewer should be \"Age: 30 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 30 Sex: F"));
      await session.step(28, "When user presses ArrowUp", () => pressKey(page, "ArrowUp"));
      await session.step(29, "Then the \"current row\" reading of markup viewer should be 1", () => readingIs(page, "current row", el("markup viewer"), 1));
      await session.step(30, "And the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Title off hides a title that is still set", async () => {
      await session.step(34, "When user sets properties of markup viewer:", () => setProperties(page, el("markup viewer"), [["Show Title","true"],["Title","Patient card"]]), [["Show Title","true"],["Title","Patient card"]]);
      await session.step(37, "Then title of markup viewer should have text \"Patient card\"", () => shouldHaveText(page, el("title of markup viewer"), "Patient card"));
      await session.step(38, "When user sets \"Show Title\" property of markup viewer to \"false\"", () => setProperty(page, "Show Title", el("markup viewer"), "false"));
      await session.step(39, "Then title of markup viewer should be hidden", () => shouldBe(page, el("title of markup viewer"), "hidden"));
      await session.step(40, "When user sets \"Show Title\" property of markup viewer to \"true\"", () => setProperty(page, "Show Title", el("markup viewer"), "true"));
      await session.step(41, "Then title of markup viewer should have text \"Patient card\"", () => shouldHaveText(page, el("title of markup viewer"), "Patient card"));
      await session.step(42, "When user sets properties of markup viewer:", () => setProperties(page, el("markup viewer"), [["Title",""],["Show Title","false"]]), [["Title",""],["Show Title","false"]]);
      await session.step(45, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit content... and OK put the new content on screen", async () => {
      await session.step(48, "Given user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(49, "Then the \"text\" reading of markup viewer should be \"Age: 26 Sex: F\"", () => readingReads(page, "text", el("markup viewer"), "Age: 26 Sex: F"));
      await session.step(50, "When user picks \"Edit content...\" from the context menu of markup viewer", () => pickFromContextMenu(page, "Edit content...", el("markup viewer")));
      await session.step(51, "Then Edit dialog should be visible", () => shouldBe(page, el("Edit dialog"), "visible"));
      await session.step(52, "When user enters \"Weight: ${WEIGHT}\" into input in Edit dialog", () => enterInto(page, "Weight: ${WEIGHT}", el("input in Edit dialog")));
      await session.step(53, "And user clicks on \"OK\" button in Edit dialog", () => clickOn(page, el("\"OK\" button in Edit dialog")));
      await session.step(54, "Then Edit dialog should be absent", () => shouldBe(page, el("Edit dialog"), "absent"));
      await session.step(55, "And \"content\" property of markup viewer should be \"Weight: ${WEIGHT}\"", () => propertyShouldBe(page, "content", el("markup viewer"), "Weight: ${WEIGHT}"));
      await session.step(56, "And the \"text\" reading of markup viewer should be \"Weight: 74.10\"", () => readingReads(page, "text", el("markup viewer"), "Weight: 74.10"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
