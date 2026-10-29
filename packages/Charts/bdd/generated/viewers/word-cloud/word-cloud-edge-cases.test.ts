/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/word-cloud/word-cloud-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.word-cloud, core#1136]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Word cloud edge cases: the word column removed", () => {
  const session = feature(test, "features/viewers/word-cloud/word-cloud-edge-cases.feature", import.meta.url);
  test("Removing the word column logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.word-cloud", "@known-failure:25", "@realizes:core#1136"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a word cloud viewer with:", () => addViewerWith(page, "word cloud", [["wordColumnName","RACE"]]), [["wordColumnName","RACE"]]);
    await session.step(16, "Then the \"words\" reading of word cloud viewer should be 4", () => readingIs(page, "words", el("word cloud viewer"), 4));
    await session.step(17, "And the \"rows shown\" reading of word cloud viewer should be 1000", () => readingIs(page, "rows shown", el("word cloud viewer"), 1000));
    await knownFailure(async () => {
      await session.step(23, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(24, "Then the table should not have a column \"RACE\"", () => hasNoColumn(page, "RACE"));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    }, 25);
  });
});
