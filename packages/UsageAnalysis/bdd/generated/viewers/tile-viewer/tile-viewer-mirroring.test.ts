/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-mirroring.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer, public#3855, public#1593]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewer} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn, renameColumn, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, hasArea, hasNoArea, hoverArea, noErrors, pickFromAreaContextMenu, pointerAway, readingAtLeast, readingDoesNotRead, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer mirrors the frame", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-mirroring.feature", import.meta.url);
  test("Tile viewer mirrors the frame", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer", "@realizes:public#3855", "@realizes:public#1593"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(18, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(19, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
    await session.step(20, "And the \"tiles\" reading of tile viewer should be at least 3", () => readingAtLeast(page, "tiles", el("tile viewer"), 3));
    await run.scenario("A field shows the grid's display string, never the raw cell, a date format picked in the grid included (public#3855)", async () => {
      await session.step(24, "Then the \"WEIGHT of row 1\" reading of tile viewer should be \"74.10\"", () => readingReads(page, "WEIGHT of row 1", el("tile viewer"), "74.10"));
      await session.step(25, "And the \"WEIGHT of row 1\" reading of tile viewer should not be \"74.1\"", () => readingDoesNotRead(page, "WEIGHT of row 1", el("tile viewer"), "74.1"));
      await session.step(26, "And the \"HEIGHT of row 1\" reading of tile viewer should be \"174.705\"", () => readingReads(page, "HEIGHT of row 1", el("tile viewer"), "174.705"));
      await session.step(27, "And the \"STARTED of row 1\" reading of tile viewer should be \"8/2/1990\"", () => readingReads(page, "STARTED of row 1", el("tile viewer"), "8/2/1990"));
      await session.step(28, "When user scrolls the grid to the \"STARTED\" column", () => scrollGridTo(page, "STARTED"));
      await session.step(29, "And user picks \"Format > dd.MM.yyyy (02.08.1990)\" from the context menu of the \"header STARTED\" area of grid", () => pickFromAreaContextMenu(page, "Format > dd.MM.yyyy (02.08.1990)", "header STARTED", el("grid")));
      await session.step(30, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(31, "Then the \"text of cell 1 of STARTED\" reading of grid should be \"02.08.1990\"", () => readingReads(page, "text of cell 1 of STARTED", el("grid"), "02.08.1990"));
      await session.step(32, "And the \"STARTED of row 1\" reading of tile viewer should be \"02.08.1990\"", () => readingReads(page, "STARTED of row 1", el("tile viewer"), "02.08.1990"));
      await session.step(33, "And the \"CONTROL of row 1\" reading of tile viewer should be \"false\"", () => readingReads(page, "CONTROL of row 1", el("tile viewer"), "false"));
      await session.step(34, "And the \"CONTROL of row 3\" reading of tile viewer should be \"true\"", () => readingReads(page, "CONTROL of row 3", el("tile viewer"), "true"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell edited in the table reaches the card", async () => {
      await session.step(38, "Then the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(39, "When user sets \"AGE\" column in row 1 to \"99\"", () => setCell(page, "AGE", 1, "99"));
      await session.step(40, "Then the value of \"AGE\" column in row 1 should be \"99\"", () => valueInRow(page, "AGE", 1, "99"));
      await session.step(41, "And the \"AGE of row 1\" reading of tile viewer should be \"99\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "99"));
      await session.step(42, "And the \"AGE of row 2\" reading of tile viewer should be \"30\"", () => readingReads(page, "AGE of row 2", el("tile viewer"), "30"));
      await session.step(43, "When user sets \"AGE\" column in row 1 to \"26\"", () => setCell(page, "AGE", 1, "26"));
      await session.step(44, "Then the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column rename carries the field, its caption and its value across", async () => {
      await session.step(48, "Then tile viewer should have a \"field AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(49, "And tile viewer should have a \"label AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "label AGE of row 1"));
      await session.step(50, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(51, "When user renames \"AGE\" column to \"AGE_YRS\"", () => renameColumn(page, "AGE", "AGE_YRS"));
      await session.step(52, "Then the table should have a column \"AGE_YRS\"", () => hasColumn(page, "AGE_YRS"));
      await session.step(53, "And the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
      await session.step(54, "And the \"fields\" reading of tile viewer should contain \"AGE_YRS\"", () => readingContains(page, "fields", el("tile viewer"), "AGE_YRS"));
      await session.step(55, "And the \"fields\" reading of tile viewer should not contain \"AGE\"", () => readingNotContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(56, "And tile viewer should have a \"field AGE_YRS of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE_YRS of row 1"));
      await session.step(57, "And tile viewer should have a \"label AGE_YRS of row 1\" area", () => hasArea(page, el("tile viewer"), "label AGE_YRS of row 1"));
      await session.step(58, "And tile viewer should not have a \"field AGE of row 1\" area", () => hasNoArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(59, "And the \"AGE_YRS of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE_YRS of row 1", el("tile viewer"), "26"));
      await session.step(60, "And the \"HEIGHT of row 1\" reading of tile viewer should be \"174.705\"", () => readingReads(page, "HEIGHT of row 1", el("tile viewer"), "174.705"));
      await session.step(61, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(62, "When user renames \"AGE_YRS\" column to \"AGE\"", () => renameColumn(page, "AGE_YRS", "AGE"));
      await session.step(63, "Then the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(64, "And the \"AGE of row 1\" reading of tile viewer should be \"26\"", () => readingReads(page, "AGE of row 1", el("tile viewer"), "26"));
      await session.step(65, "And tile viewer should have a \"field AGE of row 1\" area", () => hasArea(page, el("tile viewer"), "field AGE of row 1"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A calculated column that reaches the card shows its formatted text, not the raw number", async () => {
      await session.step(69, "When user removes \"DEMOG\" column", () => removeColumn(page, "DEMOG"));
      await session.step(70, "And user removes \"STARTED\" column", () => removeColumn(page, "STARTED"));
      await session.step(71, "Then the table should have 9 columns", () => columnCount(page, 9));
      await session.step(72, "When user adds a calculated column \"H_THIRD\" with formula \"${HEIGHT} / 3\"", () => addCalculated(page, "H_THIRD", "${HEIGHT} / 3"));
      await session.step(73, "Then the table should have 10 columns", () => columnCount(page, 10));
      await session.step(74, "And the \"fields shown\" reading of tile viewer should be 10", () => readingIs(page, "fields shown", el("tile viewer"), 10));
      await session.step(75, "And the \"fields\" reading of tile viewer should contain \"H_THIRD\"", () => readingContains(page, "fields", el("tile viewer"), "H_THIRD"));
      await session.step(76, "And tile viewer should have a \"field H_THIRD of row 1\" area", () => hasArea(page, el("tile viewer"), "field H_THIRD of row 1"));
      await session.step(77, "And the \"H_THIRD of row 1\" reading of tile viewer should be \"58.24\"", () => readingReads(page, "H_THIRD of row 1", el("tile viewer"), "58.24"));
      await session.step(78, "And the \"H_THIRD of row 1\" reading of tile viewer should not be \"58.235\"", () => readingDoesNotRead(page, "H_THIRD of row 1", el("tile viewer"), "58.235"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
      await session.step(80, "When user removes \"H_THIRD\" column", () => removeColumn(page, "H_THIRD"));
      await session.step(81, "And user adds a calculated column \"DEMOG\" with formula \"${AGE}\"", () => addCalculated(page, "DEMOG", "${AGE}"));
      await session.step(82, "And user adds a calculated column \"STARTED\" with formula \"${AGE}\"", () => addCalculated(page, "STARTED", "${AGE}"));
      await session.step(83, "Then the table should have 11 columns", () => columnCount(page, 11));
    });
    await run.scenario("A column name the card cuts off shows whole in a tooltip over its caption (public#1593)", async () => {
      await session.step(87, "Given user opens a table \"long-names\" with:", () => openTableOf(page, "long-names", [["id","measured_concentration_of_compound_in_plasma"],["1","5"],["2","6"]]), [["id","measured_concentration_of_compound_in_plasma"],["1","5"],["2","6"]]);
      await session.step(91, "And user adds a tile viewer", () => addTileViewer(page));
      await session.step(92, "When user hovers over the \"label measured_concentration_of_compound_in_plasma of row 1\" area of tile viewer", () => hoverArea(page, "label measured_concentration_of_compound_in_plasma of row 1", el("tile viewer")));
      await session.step(93, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(94, "And tooltip should contain text \"measured_concentration_of_compound_in_plasma\"", () => shouldContainText(page, el("tooltip"), "measured_concentration_of_compound_in_plasma"));
      await session.step(95, "When user hovers over the \"label id of row 1\" area of tile viewer", () => hoverArea(page, "label id of row 1", el("tile viewer")));
      await session.step(96, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(97, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
