/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaTaller, doubleClickArea, dragAreaBy, dragAreaToArea, eventFired, hasArea, hasNoArea, listenFor, noBalloons, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingDiffers, readingHigher, readingIs, readingLower, readingReads, readingSame, repainted, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid column geometry", () => {
  const session = feature(test, "features/viewers/grid/grid-columns.feature", import.meta.url);
  test("Grid column geometry", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(15, "And \"Frozen Columns\" property of grid should be \"1\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "1"));
    await session.step(16, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
    await run.scenario("Column Sizing presets order the widths Minimal, Optimal, Maximal", async () => {
      await session.step(19, "When user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 4 of AGE", el("grid")));
      await session.step(20, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(21, "When user picks \"Column Sizing > Minimal\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Minimal", "cell 4 of AGE", el("grid")));
      await session.step(22, "Then the \"column width of AGE\" reading of grid should be lower than before", () => readingLower(page, "column width of AGE", el("grid")));
      await session.step(23, "And the \"column width of RACE\" reading of grid should be lower than before", () => readingLower(page, "column width of RACE", el("grid")));
      await session.step(24, "When user picks \"Column Sizing > Maximal\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Maximal", "cell 4 of AGE", el("grid")));
      await session.step(25, "Then the \"column width of AGE\" reading of grid should be higher than before", () => readingHigher(page, "column width of AGE", el("grid")));
      await session.step(26, "And the \"column width of SEVERITY\" reading of grid should be higher than before", () => readingHigher(page, "column width of SEVERITY", el("grid")));
      await session.step(27, "When user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 4 of AGE", el("grid")));
      await session.step(28, "Then the \"column width of SEVERITY\" reading of grid should be lower than before", () => readingLower(page, "column width of SEVERITY", el("grid")));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click on the header sorts descending, then ascending, then off", async () => {
      await session.step(32, "Given user listens for \"d4-grid-rows-sorted\" event on grid", () => listenFor(page, "d4-grid-rows-sorted", el("grid")));
      await session.step(33, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(34, "Then \"d4-grid-rows-sorted\" event should have fired on grid", () => eventFired(page, "d4-grid-rows-sorted", el("grid")));
      await session.step(35, "And the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
      await session.step(36, "And the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
      await session.step(37, "And grid should have a \"cell 600 of AGE\" area", () => hasArea(page, el("grid"), "cell 600 of AGE"));
      await session.step(38, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
      await session.step(39, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(40, "Then the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(41, "And grid should have a \"cell 692 of AGE\" area", () => hasArea(page, el("grid"), "cell 692 of AGE"));
      await session.step(42, "And grid should not have a \"cell 600 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 600 of AGE"));
      await session.step(43, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(44, "Then the \"sort column\" reading of grid should be \"\"", () => readingReads(page, "sort column", el("grid"), ""));
      await session.step(45, "And grid should have a \"cell 1 of AGE\" area", () => hasArea(page, el("grid"), "cell 1 of AGE"));
      await session.step(46, "And grid should not have a \"cell 692 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 692 of AGE"));
      await session.step(47, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Sort dialog sorts by its first level", async () => {
      await session.step(51, "When user picks \"Sort...\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Sort...", "cell 4 of AGE", el("grid")));
      await session.step(52, "Then Sort Table dialog should be visible", () => shouldBe(page, el("Sort Table dialog"), "visible"));
      await session.step(53, "When user clicks on OK button in Sort Table dialog", () => clickOn(page, el("OK button in Sort Table dialog")));
      await session.step(54, "Then Sort Table dialog should be hidden", () => shouldBe(page, el("Sort Table dialog"), "hidden"));
      await session.step(55, "And the \"sort column\" reading of grid should be \"USUBJID\"", () => readingReads(page, "sort column", el("grid"), "USUBJID"));
      await session.step(56, "And the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(57, "And the value of \"AGE\" column in row 1 should be \"26\"", () => valueInRow(page, "AGE", 1, "26"));
      await session.step(58, "When user picks \"Sort > Reset\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Reset", "header AGE", el("grid")));
      await session.step(59, "Then the \"sort column\" reading of grid should be \"\"", () => readingReads(page, "sort column", el("grid"), ""));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging the header handle resizes one column and the row handle every row", async () => {
      await session.step(63, "When user drags the \"column resizer AGE\" area of grid by 40 pixels to the right", () => dragAreaBy(page, "column resizer AGE", el("grid"), 40, "right"));
      await session.step(64, "Then the \"column width of AGE\" reading of grid should be higher than before", () => readingHigher(page, "column width of AGE", el("grid")));
      await session.step(65, "And the \"column width of SEX\" reading of grid should be the same as before", () => readingSame(page, "column width of SEX", el("grid")));
      await session.step(66, "When user drags the \"row resizer 1\" area of grid by 12 pixels to the down", () => dragAreaBy(page, "row resizer 1", el("grid"), 12, "down"));
      await session.step(67, "Then the \"cell 1 of AGE\" area of grid should be taller than before", () => areaTaller(page, "cell 1 of AGE", el("grid")));
      await session.step(68, "And the \"column width of AGE\" reading of grid should be the same as before", () => readingSame(page, "column width of AGE", el("grid")));
      await session.step(69, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(70, "When user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(71, "And user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 4 of AGE", el("grid")));
      await session.step(72, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging a header onto another reorders the columns", async () => {
      await session.step(75, "When user drags the \"header HEIGHT\" area of grid to the \"header DEMOG\" area", () => dragAreaToArea(page, "header HEIGHT", el("grid"), "header DEMOG"));
      await session.step(76, "Then the \"column order\" reading of grid should differ from before", () => readingDiffers(page, "column order", el("grid")));
      await session.step(77, "And grid should have a \"header HEIGHT\" area", () => hasArea(page, el("grid"), "header HEIGHT"));
      await session.step(78, "And grid should have a \"header DEMOG\" area", () => hasArea(page, el("grid"), "header DEMOG"));
      await session.step(79, "And the table should have 11 columns", () => columnCount(page, 11));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
      await session.step(81, "When user drags the \"header HEIGHT\" area of grid to the \"header DIS_POP\" area", () => dragAreaToArea(page, "header HEIGHT", el("grid"), "header DIS_POP"));
    });
    await run.scenario("Order or Hide Columns opens on the grid and keeps its title", async () => {
      await session.step(84, "When user picks \"Order or Hide Columns...\" from the context menu of the \"cell 4 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "cell 4 of AGE", el("grid")));
      await session.step(85, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(86, "And title of Order or Hide Columns dialog should contain text \"Order or Hide Columns\"", () => shouldContainText(page, el("title of Order or Hide Columns dialog"), "Order or Hide Columns"));
      await session.step(87, "When user closes Order or Hide Columns dialog", () => close(page, el("Order or Hide Columns dialog")));
      await session.step(88, "Then Order or Hide Columns dialog should be hidden", () => shouldBe(page, el("Order or Hide Columns dialog"), "hidden"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pin Column freezes one more column and Pin Row pins two unique rows", async () => {
      await session.step(92, "Given user listens for \"d4-grid-pinned_rows-changed\" event on grid", () => listenFor(page, "d4-grid-pinned_rows-changed", el("grid")));
      await session.step(93, "When user picks \"Pin > Pin Column\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header SEX", el("grid")));
      await session.step(94, "Then \"Frozen Columns\" property of grid should be \"2\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "2"));
      await session.step(95, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 1 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 1 of USUBJID", el("grid")));
      await session.step(96, "Then \"d4-grid-pinned_rows-changed\" event should have fired on grid", () => eventFired(page, "d4-grid-pinned_rows-changed", el("grid")));
      await session.step(97, "And the \"pinned rows\" reading of grid should be 1", () => readingIs(page, "pinned rows", el("grid"), 1));
      await session.step(98, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 3 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 3 of USUBJID", el("grid")));
      await session.step(99, "Then the \"pinned rows\" reading of grid should be 2", () => readingIs(page, "pinned rows", el("grid"), 2));
      await session.step(100, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(101, "When user picks \"Pin > Unpin All Rows\" from the context menu of the \"cell 1 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Unpin All Rows", "cell 1 of USUBJID", el("grid")));
      await session.step(102, "Then the \"pinned rows\" reading of grid should be 0", () => readingIs(page, "pinned rows", el("grid"), 0));
      await session.step(103, "When user picks \"Pin > Unpin All Columns\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Unpin All Columns", "header SEX", el("grid")));
      await session.step(104, "Then \"Frozen Columns\" property of grid should be \"1\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "1"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
