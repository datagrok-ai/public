/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-core.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {pinnedPaneHidden, pinnedPaneShown} from '../../../bindings/forms.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {allOfSelected, clearSelection, filterTo, resetFilter, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, doubleClickArea, hasArea, hasNoArea, noBalloons, noErrors, pickFromAreaContextMenu, propertiesShouldBe, propertyShouldBe, readingHigher, readingIs, readingLower, readingReads, setProperty, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer core — cards, sort and pinning", () => {
  const session = feature(test, "features/viewers/forms/forms-core.feature", import.meta.url);
  test("Forms viewer core — cards, sort and pinning", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["Show Mouse Over Row","false"]]));
    await session.step(17, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(18, "And properties of forms viewer should be:", () => propertiesShouldBe(page, el("forms viewer"), [["Show Selected Rows","true"],["Show Current Row","true"],["Use Grid Sort","true"],["Renderer Size","small"],["Number Format","Same as grid"],["Color Code","true"]]));
    await session.step(25, "And the \"fields\" reading of forms viewer should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "fields", el("forms viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(26, "And the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("The current-row card shows the current row", async () => {
      await session.step(30, "When user makes row 13 current", () => makeRowCurrent(page, 13));
      await session.step(31, "Then the \"current record\" reading of forms viewer should be 13", () => readingIs(page, "current record", el("forms viewer"), 13));
      await session.step(32, "And the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000900008\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000900008"));
      await session.step(33, "And the \"AGE of current card\" reading of forms viewer should be \"43\"", () => readingReads(page, "AGE of current card", el("forms viewer"), "43"));
      await session.step(34, "And the \"cards\" reading of forms viewer should be 1", () => readingIs(page, "cards", el("forms viewer"), 1));
      await session.step(35, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(36, "Then the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000300003"));
      await session.step(37, "And the \"AGE of current card\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of current card", el("forms viewer"), "26"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("One card per selected row, in table order", async () => {
      await session.step(41, "When user selects rows where \"USUBJID\" is one of \"X0273T21000400001, X0273T21000500008, X0273T21001500015\"", () => selectWhereOneOf(page, "USUBJID", "X0273T21000400001, X0273T21000500008, X0273T21001500015"));
      await session.step(42, "Then the \"cards\" reading of forms viewer should be 4", () => readingIs(page, "cards", el("forms viewer"), 4));
      await session.step(43, "And the \"records shown\" reading of forms viewer should be 4", () => readingIs(page, "records shown", el("forms viewer"), 4));
      await session.step(44, "And the \"USUBJID of card 2\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of card 2", el("forms viewer"), "X0273T21000400001"));
      await session.step(45, "And the \"USUBJID of card 3\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of card 3", el("forms viewer"), "X0273T21000500008"));
      await session.step(46, "And the \"USUBJID of card 4\" reading of forms viewer should be \"X0273T21001500015\"", () => readingReads(page, "USUBJID of card 4", el("forms viewer"), "X0273T21001500015"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filtered-out selected row loses its card", async () => {
      await session.step(50, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(51, "Then the \"cards\" reading of forms viewer should be 2", () => readingIs(page, "cards", el("forms viewer"), 2));
      await session.step(52, "And the \"USUBJID of card 2\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of card 2", el("forms viewer"), "X0273T21000500008"));
      await session.step(53, "When user resets the filter", () => resetFilter(page));
      await session.step(54, "Then the \"cards\" reading of forms viewer should be 4", () => readingIs(page, "cards", el("forms viewer"), 4));
      await session.step(55, "And the \"USUBJID of card 2\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of card 2", el("forms viewer"), "X0273T21000400001"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click on a label advances the same sort cycle", async () => {
      await session.step(59, "When user sets \"Sort By\" property of forms viewer to \"AGE\"", () => setProperty(page, "Sort By", el("forms viewer"), "AGE"));
      await session.step(60, "Then the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
      await session.step(61, "And the \"sort direction\" reading of forms viewer should be \"↓\"", () => readingReads(page, "sort direction", el("forms viewer"), "↓"));
      await session.step(62, "When user double-clicks on the \"label AGE\" area of forms viewer", () => doubleClickArea(page, "label AGE", el("forms viewer")));
      await session.step(63, "Then the \"sort direction\" reading of forms viewer should be \"↑\"", () => readingReads(page, "sort direction", el("forms viewer"), "↑"));
      await session.step(64, "And the \"sort column\" reading of forms viewer should be \"AGE\"", () => readingReads(page, "sort column", el("forms viewer"), "AGE"));
      await session.step(65, "When user double-clicks on the \"label HEIGHT\" area of forms viewer", () => doubleClickArea(page, "label HEIGHT", el("forms viewer")));
      await session.step(66, "Then the \"sort column\" reading of forms viewer should be \"\"", () => readingReads(page, "sort column", el("forms viewer"), ""));
      await session.step(67, "And forms viewer should not have a \"sort indicator HEIGHT\" area", () => hasNoArea(page, el("forms viewer"), "sort indicator HEIGHT"));
      await session.step(68, "When user double-clicks on the \"label HEIGHT\" area of forms viewer", () => doubleClickArea(page, "label HEIGHT", el("forms viewer")));
      await session.step(69, "Then the \"sort column\" reading of forms viewer should be \"HEIGHT\"", () => readingReads(page, "sort column", el("forms viewer"), "HEIGHT"));
      await session.step(70, "And the \"sort direction\" reading of forms viewer should be \"↓\"", () => readingReads(page, "sort direction", el("forms viewer"), "↓"));
      await session.step(71, "When user sets \"Sort By\" property of forms viewer to \"\"", () => setProperty(page, "Sort By", el("forms viewer"), ""));
      await session.step(72, "Then the \"sort column\" reading of forms viewer should be \"\"", () => readingReads(page, "sort column", el("forms viewer"), ""));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid's sort orders the cards and marks the label", async () => {
      await session.step(76, "When user picks \"Sort > Ascending\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Ascending", "header HEIGHT", el("grid")));
      await session.step(77, "Then the \"sort column\" reading of grid should be \"HEIGHT\"", () => readingReads(page, "sort column", el("grid"), "HEIGHT"));
      await session.step(78, "And the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(79, "And the \"sort column\" reading of forms viewer should be \"HEIGHT\"", () => readingReads(page, "sort column", el("forms viewer"), "HEIGHT"));
      await session.step(80, "And the \"sort direction\" reading of forms viewer should be \"↑\"", () => readingReads(page, "sort direction", el("forms viewer"), "↑"));
      await session.step(81, "And forms viewer should have a \"sort indicator HEIGHT\" area", () => hasArea(page, el("forms viewer"), "sort indicator HEIGHT"));
      await session.step(82, "And the \"USUBJID of card 2\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of card 2", el("forms viewer"), "X0273T21000500008"));
      await session.step(83, "And the \"USUBJID of card 3\" reading of forms viewer should be \"X0273T21001500015\"", () => readingReads(page, "USUBJID of card 3", el("forms viewer"), "X0273T21001500015"));
      await session.step(84, "And the \"USUBJID of card 4\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of card 4", el("forms viewer"), "X0273T21000400001"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sort By orders the cards without touching the grid's own sort", async () => {
      await session.step(88, "When user sets \"Sort By\" property of forms viewer to \"WEIGHT\"", () => setProperty(page, "Sort By", el("forms viewer"), "WEIGHT"));
      await session.step(89, "Then the \"sort column\" reading of forms viewer should be \"WEIGHT\"", () => readingReads(page, "sort column", el("forms viewer"), "WEIGHT"));
      await session.step(90, "And the \"sort direction\" reading of forms viewer should be \"↓\"", () => readingReads(page, "sort direction", el("forms viewer"), "↓"));
      await session.step(91, "And forms viewer should have a \"sort indicator WEIGHT\" area", () => hasArea(page, el("forms viewer"), "sort indicator WEIGHT"));
      await session.step(92, "And forms viewer should not have a \"sort indicator HEIGHT\" area", () => hasNoArea(page, el("forms viewer"), "sort indicator HEIGHT"));
      await session.step(93, "And the \"sort column\" reading of grid should be \"HEIGHT\"", () => readingReads(page, "sort column", el("grid"), "HEIGHT"));
      await session.step(94, "And the \"USUBJID of card 2\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of card 2", el("forms viewer"), "X0273T21000400001"));
      await session.step(95, "And the \"USUBJID of card 3\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of card 3", el("forms viewer"), "X0273T21000500008"));
      await session.step(96, "And the \"USUBJID of card 4\" reading of forms viewer should be \"X0273T21001500015\"", () => readingReads(page, "USUBJID of card 4", el("forms viewer"), "X0273T21001500015"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Use Grid Sort off stops the cards from mirroring the grid", async () => {
      await session.step(100, "When user sets \"Sort By\" property of forms viewer to \"\"", () => setProperty(page, "Sort By", el("forms viewer"), ""));
      await session.step(101, "Then the \"sort column\" reading of forms viewer should be \"HEIGHT\"", () => readingReads(page, "sort column", el("forms viewer"), "HEIGHT"));
      await session.step(102, "When user sets \"Use Grid Sort\" property of forms viewer to \"false\"", () => setProperty(page, "Use Grid Sort", el("forms viewer"), "false"));
      await session.step(103, "Then the \"sort column\" reading of forms viewer should be \"\"", () => readingReads(page, "sort column", el("forms viewer"), ""));
      await session.step(104, "And the \"USUBJID of card 2\" reading of forms viewer should be \"X0273T21000400001\"", () => readingReads(page, "USUBJID of card 2", el("forms viewer"), "X0273T21000400001"));
      await session.step(105, "When user sets \"Use Grid Sort\" property of forms viewer to \"true\"", () => setProperty(page, "Use Grid Sort", el("forms viewer"), "true"));
      await session.step(106, "Then the \"sort column\" reading of forms viewer should be \"HEIGHT\"", () => readingReads(page, "sort column", el("forms viewer"), "HEIGHT"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pin Row moves a card into the pinned pane and Unpin Row brings it back", async () => {
      await session.step(110, "When user picks \"Pin Row\" from the context menu of the \"field USUBJID of card 2\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field USUBJID of card 2", el("forms viewer")));
      await session.step(111, "Then the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
      await session.step(112, "And the pinned pane of forms viewer should be shown", () => pinnedPaneShown(page));
      await session.step(113, "And the \"USUBJID of pinned card 1\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of pinned card 1", el("forms viewer"), "X0273T21000500008"));
      await session.step(114, "And the \"cards\" reading of forms viewer should be lower than before", () => readingLower(page, "cards", el("forms viewer")));
      await session.step(115, "And all rows where \"USUBJID\" is \"X0273T21000500008\" should be selected", () => allOfSelected(page, "USUBJID", "X0273T21000500008"));
      await session.step(116, "And \"pinnedRowValues\" property of forms viewer should be \"X0273T21000500008\"", () => propertyShouldBe(page, "pinnedRowValues", el("forms viewer"), "X0273T21000500008"));
      await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(118, "When user picks \"Unpin Row\" from the context menu of the \"pinned card 1\" area of forms viewer", () => pickFromAreaContextMenu(page, "Unpin Row", "pinned card 1", el("forms viewer")));
      await session.step(119, "Then the \"pinned records\" reading of forms viewer should be 0", () => readingIs(page, "pinned records", el("forms viewer"), 0));
      await session.step(120, "And the pinned pane of forms viewer should be hidden", () => pinnedPaneHidden(page));
      await session.step(121, "And the \"cards\" reading of forms viewer should be higher than before", () => readingHigher(page, "cards", el("forms viewer")));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pinning through a non-unique value warns and pins all the same", async () => {
      await session.step(125, "When user picks \"Pin Row\" from the context menu of the \"field SEX of card 2\" area of forms viewer", () => pickFromAreaContextMenu(page, "Pin Row", "field SEX of card 2", el("forms viewer")));
      await session.step(126, "Then a warning balloon containing \"non-unique value\" should have been shown", () => warningBalloonText(page, "non-unique value"));
      await session.step(127, "And the \"pinned records\" reading of forms viewer should be 1", () => readingIs(page, "pinned records", el("forms viewer"), 1));
      await session.step(128, "And \"pinnedRowValues\" property of forms viewer should be \"M\"", () => propertyShouldBe(page, "pinnedRowValues", el("forms viewer"), "M"));
      await session.step(129, "When user picks \"Unpin Row\" from the context menu of the \"pinned card 1\" area of forms viewer", () => pickFromAreaContextMenu(page, "Unpin Row", "pinned card 1", el("forms viewer")));
      await session.step(130, "Then the \"pinned records\" reading of forms viewer should be 0", () => readingIs(page, "pinned records", el("forms viewer"), 0));
      await session.step(131, "And the pinned pane of forms viewer should be hidden", () => pinnedPaneHidden(page));
      await session.step(132, "And the \"cards\" reading of forms viewer should be 4", () => readingIs(page, "cards", el("forms viewer"), 4));
      await session.step(133, "When user clears the row selection", () => clearSelection(page));
      await session.step(134, "Then the \"cards\" reading of forms viewer should be 1", () => readingIs(page, "cards", el("forms viewer"), 1));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
