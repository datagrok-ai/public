/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-mirroring.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addTileViewer} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn, renameColumn, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hasNoArea, noErrors, readingAtLeast, readingDoesNotRead, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer mirrors the frame", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-mirroring.feature", import.meta.url);
  test("Tile viewer mirrors the frame", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(15, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(16, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
    await session.step(17, "And the \"tiles\" reading of tile viewer should be at least 3", () => readingAtLeast(page, "tiles", el("tile viewer"), 3));
    await run.scenario("A field shows the grid's display string, never the raw cell", async () => {
      await session.step(20, "Then the \"WEIGHT of row 1\" reading of tile viewer should be \"74.10\"", () => readingReads(page, "WEIGHT of row 1", el("tile viewer"), "74.10"));
      await session.step(21, "And the \"WEIGHT of row 1\" reading of tile viewer should not be \"74.1\"", () => readingDoesNotRead(page, "WEIGHT of row 1", el("tile viewer"), "74.1"));
      await session.step(22, "And the \"HEIGHT of row 1\" reading of tile viewer should be \"174.705\"", () => readingReads(page, "HEIGHT of row 1", el("tile viewer"), "174.705"));
      await session.step(23, "And the \"STARTED of row 1\" reading of tile viewer should be \"8/2/1990\"", () => readingReads(page, "STARTED of row 1", el("tile viewer"), "8/2/1990"));
      await session.step(24, "And the \"CONTROL of row 1\" reading of tile viewer should be \"false\"", () => readingReads(page, "CONTROL of row 1", el("tile viewer"), "false"));
      await session.step(25, "And the \"CONTROL of row 3\" reading of tile viewer should be \"true\"", () => readingReads(page, "CONTROL of row 3", el("tile viewer"), "true"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell edited in the table reaches the card", async () => {
      await session.step(29, "Then the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(30, "When user sets \"AGE\" column in row 1 to \"99\"", () => setCell(page, "AGE", 1, "99"));
      await session.step(31, "Then the value of \"AGE\" column in row 1 should be \"99\"", () => valueInRow(page, "AGE", 1, "99"));
      await session.step(32, "And the \"AGE of row 1\" reading of tile viewer should be \"99\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "99"));
      await session.step(33, "And the \"AGE of row 2\" reading of tile viewer should be \"30\"", () => readingReads(page, "AGE of row 2", el("tile viewer"), "30"));
      await session.step(34, "When user sets \"AGE\" column in row 1 to \"26\"", () => setCell(page, "AGE", 1, "26"));
      await session.step(35, "Then the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column rename carries the field, its caption and its value across", async () => {
      await session.step(39, "Then tile viewer should have a \"field AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(40, "And tile viewer should have a \"label AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "label AGE of row 1"));
      await session.step(41, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(42, "When user renames \"AGE\" column to \"AGE_YRS\"", () => renameColumn(page, "AGE", "AGE_YRS"));
      await session.step(43, "Then the table should have a column \"AGE_YRS\"", () => hasColumn(page, "AGE_YRS"));
      await session.step(44, "And the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
      await session.step(45, "And the \"fields\" reading of tile viewer should contain \"AGE_YRS\"", () => readingContains(page, "fields", el("tile viewer"), "AGE_YRS"));
      await session.step(46, "And the \"fields\" reading of tile viewer should not contain \"AGE\"", () => readingNotContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(47, "And tile viewer should have a \"field AGE_YRS of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE_YRS of row 1"));
      await session.step(48, "And tile viewer should have a \"label AGE_YRS of row 1\" area", () => hasArea(page, el("tile viewer"), "label AGE_YRS of row 1"));
      await session.step(49, "And tile viewer should not have a \"field AGE of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(50, "And the \"AGE_YRS of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE_YRS of row 1", el("tile viewer"), "26"));
      await session.step(51, "And the \"HEIGHT of row 1\" reading of tile viewer should be \"174.705\"", () => readingReads(page, "HEIGHT of row 1", el("tile viewer"), "174.705"));
      await session.step(52, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(53, "When user renames \"AGE_YRS\" column to \"AGE\"", () => renameColumn(page, "AGE_YRS", "AGE"));
      await session.step(54, "Then the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(55, "And the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(56, "And tile viewer should have a \"field AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A calculated column that reaches the card shows its formatted text, not the raw number", async () => {
      await session.step(60, "When user removes \"DEMOG\" column", () => removeColumn(page, "DEMOG"));
      await session.step(61, "And user removes \"STARTED\" column", () => removeColumn(page, "STARTED"));
      await session.step(62, "Then the table should have 9 columns", () => columnCount(page, 9));
      await session.step(63, "When user adds a calculated column \"H_THIRD\" with formula \"${HEIGHT} / 3\"", () => addCalculated(page, "H_THIRD", "${HEIGHT} / 3"));
      await session.step(64, "Then the table should have 10 columns", () => columnCount(page, 10));
      await session.step(65, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(66, "And the \"fields\" reading of tile viewer should contain \"H_THIRD\"", () => readingContains(page, "fields", el("tile viewer"), "H_THIRD"));
      await session.step(67, "And tile viewer should have a \"field H_THIRD of row 1\" area", () => hasArea(page, el("tile viewer"), "field H_THIRD of row 1"));
      await session.step(68, "And the \"H_THIRD of row 1\" reading of tile viewer should be \"58.24\"", () => readingReads(page, "H_THIRD of row 1", el("tile viewer"), "58.24"));
      await session.step(69, "And the \"H_THIRD of row 1\" reading of tile viewer should not be \"58.235\"", () => readingDoesNotRead(page, "H_THIRD of row 1", el("tile viewer"), "58.235"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
      await session.step(71, "When user removes \"H_THIRD\" column", () => removeColumn(page, "H_THIRD"));
      await session.step(72, "And user adds a calculated column \"DEMOG\" with formula \"${AGE}\"", () => addCalculated(page, "DEMOG", "${AGE}"));
      await session.step(73, "And user adds a calculated column \"STARTED\" with formula \"${AGE}\"", () => addCalculated(page, "STARTED", "${AGE}"));
      await session.step(74, "Then the table should have 11 columns", () => columnCount(page, 11));
    });
    run.finish();
  });
});
