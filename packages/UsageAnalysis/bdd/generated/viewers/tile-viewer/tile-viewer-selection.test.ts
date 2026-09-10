/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addTileViewer} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, noneOfSelected, noneSelected, onlyOfAnySelected, onlyOfSelected, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, noErrors, propertyShouldBe, readingAtLeast, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {everyTileShows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer current row and selection", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-selection.feature", import.meta.url);
  test("Tile viewer current row and selection", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(18, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(19, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(20, "And the \"tiles\" reading of tile viewer should be at least 5", () => readingAtLeast(page, "tiles", el("tile viewer"), 5));
    await run.scenario("Nothing is selected and the first row is current", async () => {
      await session.step(23, "Then the \"current row\" reading of tile viewer should be 1", () => readingIs(page, "current row", el("tile viewer"), 1));
      await session.step(24, "And the \"rows selected\" reading of tile viewer should be 0", () => readingIs(page, "rows selected", el("tile viewer"), 0));
      await session.step(25, "And no rows should be selected", () => noneSelected(page));
      await session.step(26, "And the \"selected rows shown\" reading of tile viewer should be \"true\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "true"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A plain click makes the card's row current and selects nothing", async () => {
      await session.step(30, "When user clicks on the \"tile of row 4\" area of tile viewer", () => clickArea(page, "tile of row 4", el("tile viewer")));
      await session.step(31, "Then the \"current row\" reading of tile viewer should be 4", () => readingIs(page, "current row", el("tile viewer"), 4));
      await session.step(32, "And row 4 should be current", () => currentRowIs(page, 4));
      await session.step(33, "And the \"rows selected\" reading of tile viewer should be 0", () => readingIs(page, "rows selected", el("tile viewer"), 0));
      await session.step(34, "And no rows should be selected", () => noneSelected(page));
      await session.step(35, "When user clicks on the \"tile of row 2\" area of tile viewer", () => clickArea(page, "tile of row 2", el("tile viewer")));
      await session.step(36, "Then the \"current row\" reading of tile viewer should be 2", () => readingIs(page, "current row", el("tile viewer"), 2));
      await session.step(37, "And row 2 should be current", () => currentRowIs(page, 2));
      await session.step(38, "And no rows should be selected", () => noneSelected(page));
      await session.step(39, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(40, "Then the \"current row\" reading of tile viewer should be 1", () => readingIs(page, "current row", el("tile viewer"), 1));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control-click adds the card's row to the selection and takes it back out", async () => {
      await session.step(44, "Given no rows should be selected", () => noneSelected(page));
      await session.step(45, "When user clicks on the \"tile of row 5\" area of tile viewer holding Control", () => clickAreaHolding(page, "tile of row 5", el("tile viewer"), "Control"));
      await session.step(46, "Then the \"rows selected\" reading of tile viewer should be 1", () => readingIs(page, "rows selected", el("tile viewer"), 1));
      await session.step(47, "And only rows where \"USUBJID\" is \"X0273T21000500006\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000500006"));
      await session.step(48, "When user clicks on the \"tile of row 5\" area of tile viewer holding Control", () => clickAreaHolding(page, "tile of row 5", el("tile viewer"), "Control"));
      await session.step(49, "Then the \"rows selected\" reading of tile viewer should be 0", () => readingIs(page, "rows selected", el("tile viewer"), 0));
      await session.step(50, "And no rows should be selected", () => noneSelected(page));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift-click adds one row, not the range up to it", async () => {
      await session.step(54, "Given no rows should be selected", () => noneSelected(page));
      await session.step(55, "When user clicks on the \"tile of row 1\" area of tile viewer", () => clickArea(page, "tile of row 1", el("tile viewer")));
      await session.step(56, "Then the \"current row\" reading of tile viewer should be 1", () => readingIs(page, "current row", el("tile viewer"), 1));
      await session.step(57, "When user clicks on the \"tile of row 3\" area of tile viewer holding Shift", () => clickAreaHolding(page, "tile of row 3", el("tile viewer"), "Shift"));
      await session.step(58, "Then the \"rows selected\" reading of tile viewer should be 1", () => readingIs(page, "rows selected", el("tile viewer"), 1));
      await session.step(59, "And only rows where \"USUBJID\" is \"X0273T21000400001\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000400001"));
      await session.step(60, "And no rows where \"USUBJID\" is \"X0273T21000300005\" should be selected", () => noneOfSelected(page, "USUBJID", "X0273T21000300005"));
      await session.step(61, "When user clears the row selection", () => clearSelection(page));
      await session.step(62, "Then no rows should be selected", () => noneSelected(page));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control+Shift-click takes one row out of the selection and leaves the rest", async () => {
      await session.step(66, "Given no rows should be selected", () => noneSelected(page));
      await session.step(67, "When user clicks on the \"tile of row 5\" area of tile viewer holding Control", () => clickAreaHolding(page, "tile of row 5", el("tile viewer"), "Control"));
      await session.step(68, "And user clicks on the \"tile of row 6\" area of tile viewer holding Control", () => clickAreaHolding(page, "tile of row 6", el("tile viewer"), "Control"));
      await session.step(69, "Then the \"rows selected\" reading of tile viewer should be 2", () => readingIs(page, "rows selected", el("tile viewer"), 2));
      await session.step(70, "And only rows where \"USUBJID\" is one of \"X0273T21000500006, X0273T21000500008\" should be selected", () => onlyOfAnySelected(page, "USUBJID", "X0273T21000500006, X0273T21000500008"));
      await session.step(71, "When user clicks on the \"tile of row 5\" area of tile viewer holding Control+Shift", () => clickAreaHolding(page, "tile of row 5", el("tile viewer"), "Control+Shift"));
      await session.step(72, "Then the \"rows selected\" reading of tile viewer should be 1", () => readingIs(page, "rows selected", el("tile viewer"), 1));
      await session.step(73, "And only rows where \"USUBJID\" is \"X0273T21000500008\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000500008"));
      await session.step(74, "And no rows where \"USUBJID\" is \"X0273T21000500006\" should be selected", () => noneOfSelected(page, "USUBJID", "X0273T21000500006"));
      await session.step(75, "When user clears the row selection", () => clearSelection(page));
      await session.step(76, "Then no rows should be selected", () => noneSelected(page));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Selected Rows off neutralises the highlight and keeps the selection", async () => {
      await session.step(80, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(81, "Then the \"rows selected\" reading of tile viewer should be 15", () => readingIs(page, "rows selected", el("tile viewer"), 15));
      await session.step(82, "And the \"selected rows shown\" reading of tile viewer should be \"true\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "true"));
      await session.step(83, "When user sets \"Show Selected Rows\" property of tile viewer to \"false\"", () => setProperty(page, "Show Selected Rows", el("tile viewer"), "false"));
      await session.step(84, "Then the \"selected rows shown\" reading of tile viewer should be \"false\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "false"));
      await session.step(85, "And the \"rows selected\" reading of tile viewer should be 15", () => readingIs(page, "rows selected", el("tile viewer"), 15));
      await session.step(86, "And 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(87, "When user sets \"Show Selected Rows\" property of tile viewer to \"true\"", () => setProperty(page, "Show Selected Rows", el("tile viewer"), "true"));
      await session.step(88, "Then the \"selected rows shown\" reading of tile viewer should be \"true\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "true"));
      await session.step(89, "And 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source = Selected suppresses the highlight whatever the property says", async () => {
      await session.step(93, "Then \"Show Selected Rows\" property of tile viewer should be \"true\"", () => propertyShouldBe(page, "Show Selected Rows", el("tile viewer"), "true"));
      await session.step(94, "And the \"selected rows shown\" reading of tile viewer should be \"true\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "true"));
      await session.step(95, "When user sets \"Row Source\" property of tile viewer to \"Selected\"", () => setProperty(page, "Row Source", el("tile viewer"), "Selected"));
      await session.step(96, "Then tile viewer should show 15 rows", () => showsRows(page, el("tile viewer"), 15));
      await session.step(97, "And \"Show Selected Rows\" property of tile viewer should be \"true\"", () => propertyShouldBe(page, "Show Selected Rows", el("tile viewer"), "true"));
      await session.step(98, "And the \"selected rows shown\" reading of tile viewer should be \"false\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "false"));
      await session.step(99, "And every tile of tile viewer should show \"Asian\" in \"RACE\"", () => everyTileShows(page, el("tile viewer"), "Asian", "RACE"));
      await session.step(100, "When user sets \"Row Source\" property of tile viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("tile viewer"), "Filtered"));
      await session.step(101, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
      await session.step(102, "And the \"selected rows shown\" reading of tile viewer should be \"true\"", () => readingReads(page, "selected rows shown", el("tile viewer"), "true"));
      await session.step(103, "When user clears the row selection", () => clearSelection(page));
      await session.step(104, "Then no rows should be selected", () => noneSelected(page));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
