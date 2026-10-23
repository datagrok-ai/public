/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-stability.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend]
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterTo, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragAcrossArea, legendDocked, legendInCorner, legendLists, legendMiniIcon, legendPlacedAsBefore, legendSameItems, legendSlot, noErrors, pickFromContextMenu, pointerAway, readingDiffers, readingLower, resizeTo, restoreSize, setProperty, showsRows, takeSnapshot, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A legend that stays where it is put", () => {
  const session = feature(test, "features/viewers/legend/legend-stability.feature", import.meta.url);
  test("A legend that stays where it is put", {tag: ["@journey", "@viewers", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Zoom and Filter","no action"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Zoom and Filter","no action"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(30, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
    await session.step(31, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await run.scenario("A wheel zoom and a pan leave the docked legend where it was", async () => {
      await session.step(34, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(35, "And user scrolls the mouse wheel up 3 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "up", 3, "view", el("scatter plot viewer")));
      await session.step(36, "Then the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(37, "And the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(38, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(39, "And the legend of scatter plot viewer should be in the \"right\" slot", () => legendSlot(page, el("scatter plot viewer"), "right"));
      await session.step(40, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(41, "And user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(42, "Then the \"x axis min\" reading of scatter plot viewer should differ from before", () => readingDiffers(page, "x axis min", el("scatter plot viewer")));
      await session.step(43, "And the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(44, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(45, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(46, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A wheel zoom and a pan leave a corner legend in its corner", async () => {
      await session.step(49, "When user sets \"Legend Position\" property of scatter plot viewer to \"RightBottom\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "RightBottom"));
      await session.step(50, "Then the legend of scatter plot viewer should be in a corner", () => legendInCorner(page, el("scatter plot viewer")));
      await session.step(51, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(52, "And user scrolls the mouse wheel up 3 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "up", 3, "view", el("scatter plot viewer")));
      await session.step(53, "Then the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(54, "And the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(55, "And the legend of scatter plot viewer should be in the \"rightBottom\" slot", () => legendSlot(page, el("scatter plot viewer"), "rightBottom"));
      await session.step(56, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(57, "And user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(58, "Then the \"x axis min\" reading of scatter plot viewer should differ from before", () => readingDiffers(page, "x axis min", el("scatter plot viewer")));
      await session.step(59, "And the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(60, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(61, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(62, "And user sets \"Legend Position\" property of scatter plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Right"));
      await session.step(63, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A corner legend folded by its chevron stays folded through Position Auto and resizes", async () => {
      await session.step(67, "When user sets \"Legend Position\" property of scatter plot viewer to \"RightBottom\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "RightBottom"));
      await session.step(68, "Then the legend of scatter plot viewer should be in a corner", () => legendInCorner(page, el("scatter plot viewer")));
      await session.step(69, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(70, "And user hovers over legend of scatter plot viewer", () => hoverOver(page, el("legend of scatter plot viewer")));
      await session.step(71, "And user clicks on legend close chevron of scatter plot viewer", () => clickOn(page, el("legend close chevron of scatter plot viewer")));
      await session.step(72, "Then the legend of scatter plot viewer should be collapsed to the mini icon", () => legendMiniIcon(page, el("scatter plot viewer")));
      await session.step(73, "When user sets \"Legend Position\" property of scatter plot viewer to \"Auto\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Auto"));
      await session.step(74, "And user resizes scatter plot viewer to 380 by 300", () => resizeTo(page, el("scatter plot viewer"), 380, 300));
      await session.step(75, "Then the legend of scatter plot viewer should be collapsed to the mini icon", () => legendMiniIcon(page, el("scatter plot viewer")));
      await session.step(76, "When user resizes scatter plot viewer to 950 by 700", () => resizeTo(page, el("scatter plot viewer"), 950, 700));
      await session.step(77, "Then the legend of scatter plot viewer should be collapsed to the mini icon", () => legendMiniIcon(page, el("scatter plot viewer")));
      await session.step(78, "When user resizes scatter plot viewer to 700 by 620", () => resizeTo(page, el("scatter plot viewer"), 700, 620));
      await session.step(79, "Then the legend of scatter plot viewer should be collapsed to the mini icon", () => legendMiniIcon(page, el("scatter plot viewer")));
      await session.step(80, "And mini legend icon of scatter plot viewer should be visible", () => shouldBe(page, el("mini legend icon of scatter plot viewer"), "visible"));
      await session.step(81, "And legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
      await session.step(82, "When user clicks on mini legend icon of scatter plot viewer", () => clickOn(page, el("mini legend icon of scatter plot viewer")));
      await session.step(83, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
      await session.step(84, "And mini legend icon of scatter plot viewer should be hidden", () => shouldBe(page, el("mini legend icon of scatter plot viewer"), "hidden"));
      await session.step(85, "When user restores the size of scatter plot viewer", () => restoreSize(page, el("scatter plot viewer")));
      await session.step(86, "And user sets \"Legend Position\" property of scatter plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Right"));
      await session.step(87, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A legend column of a thousand categories docks on the right and lists the rows drawn", async () => {
      await session.step(91, "When user sets \"Color\" property of scatter plot viewer to \"USUBJID\"", () => setProperty(page, "Color", el("scatter plot viewer"), "USUBJID"));
      await session.step(92, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(93, "And the legend of scatter plot viewer should be in the \"right\" slot", () => legendSlot(page, el("scatter plot viewer"), "right"));
      await session.step(94, "And the legend of scatter plot viewer should list 872 items", () => legendLists(page, el("scatter plot viewer"), 872));
      await session.step(95, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(96, "When user sets \"Legend Position\" property of scatter plot viewer to \"Auto\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Auto"));
      await session.step(97, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(98, "And the legend of scatter plot viewer should list 872 items", () => legendLists(page, el("scatter plot viewer"), 872));
      await session.step(99, "When user sets \"Legend Position\" property of scatter plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Right"));
      await session.step(100, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(101, "Then the legend of scatter plot viewer should list 499 items", () => legendLists(page, el("scatter plot viewer"), 499));
      await session.step(102, "And scatter plot viewer should show 499 rows", () => showsRows(page, el("scatter plot viewer"), 499));
      await session.step(103, "When user resets the filter", () => resetFilter(page));
      await session.step(104, "And user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(105, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
