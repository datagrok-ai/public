/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/word-cloud/word-cloud-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.word-cloud, core#1136]
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
import {clickOn, enterInto, finishedUpdating, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPasses, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Word cloud edge cases: no rows, Rotation Step below its minimum, the word column removed", () => {
  const session = feature(test, "features/viewers/word-cloud/word-cloud-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves no row counted and no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.word-cloud"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"]]), [["wordColumnName","RACE"]]);
    await session.step(19, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    await session.step(20, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
    await session.step(23, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(24, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(25, "And the \"rows shown\" reading of word cloud viewer should be 0", () => readingIs(page, "rows shown", el("word cloud viewer"), 0));
    await session.step(26, "And word cloud viewer should have finished updating", () => finishedUpdating(page, el("word cloud viewer")));
    await session.step(27, "When user resets the filter", () => resetFilter(page));
    await session.step(28, "Then the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
    await session.step(29, "And the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    await session.step(30, "And the \"rows of word \\\"Caucasian\\\"\" reading of word cloud viewer should be 896", () => readingIs(page, "rows of word \"Caucasian\"", el("word cloud viewer"), 896));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("Rotation Step 0 typed in the context panel keeps the words drawn without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.word-cloud", "@known-failure:43", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"]]), [["wordColumnName","RACE"]]);
    await session.step(19, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    await session.step(20, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
    await knownFailure(async () => {
      await session.step(38, "When user clicks on settings icon of word cloud viewer", () => clickOn(page, el("settings icon of word cloud viewer")));
      await session.step(39, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(40, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(41, "When user enters \"0\" into \"Rotation Step\" property in context panel", () => enterInto(page, "0", el("\"Rotation Step\" property in context panel")));
      await session.step(42, "Then word cloud viewer should have finished updating", () => finishedUpdating(page, el("word cloud viewer")));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
      await session.step(44, "And the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    }, 43);
  });
  test("Removing the word column logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.word-cloud", "@known-failure:52", "@realizes:core#1136"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"]]), [["wordColumnName","RACE"]]);
    await session.step(19, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    await session.step(20, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
    await knownFailure(async () => {
      await session.step(50, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(51, "Then the table should not have a column \"RACE\"", () => hasNoColumn(page, "RACE"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    }, 52);
  });
});
