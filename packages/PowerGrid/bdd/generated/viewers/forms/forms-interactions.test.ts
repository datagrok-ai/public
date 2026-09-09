/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-interactions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {chordClick, currentColumnIs} from '../../../bindings/forms.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {allOfSelected, currentRowValue, noneOfSelected, selectWhereOneOf, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, propertiesShouldBe, readingHigher, readingIs, readingLower, readingReads, readingSame, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer interactions and row binding", () => {
  const session = feature(test, "features/viewers/forms/forms-interactions.feature", import.meta.url);
  test("Forms viewer interactions and row binding", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(14, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(15, "And properties of forms viewer should be:", () => propertiesShouldBe(page, el("forms viewer"), [["Show Current Row","true"],["Show Mouse Over Row","true"],["Show Selected Rows","true"]]));
    await session.step(19, "When user makes row 6 current", () => makeRowCurrent(page, 6));
    await session.step(20, "And user selects rows where \"USUBJID\" is one of \"X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015\"", () => selectWhereOneOf(page, "USUBJID", "X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015"));
    await session.step(21, "Then the \"cards\" reading of forms viewer should be 6", () => readingIs(page, "cards", el("forms viewer"), 6));
    await session.step(22, "And the \"records shown\" reading of forms viewer should be 5", () => readingIs(page, "records shown", el("forms viewer"), 5));
    await session.step(23, "And the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000500008"));
    await run.scenario("A click on a card makes its row current", async () => {
      await session.step(26, "When user clicks on the \"card 5\" area of forms viewer", () => clickArea(page, "card 5", el("forms viewer")));
      await session.step(27, "Then \"USUBJID\" of the current row should be \"X0273T21000900003\"", () => currentRowValue(page, "USUBJID", "X0273T21000900003"));
      await session.step(28, "And the \"current record\" reading of forms viewer should be 12", () => readingIs(page, "current record", el("forms viewer"), 12));
      await session.step(29, "And the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000900003\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000900003"));
      await session.step(30, "When user makes row 6 current", () => makeRowCurrent(page, 6));
      await session.step(31, "Then the \"current record\" reading of forms viewer should be 6", () => readingIs(page, "current record", el("forms viewer"), 6));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a field sets the current cell", async () => {
      await session.step(35, "When user clicks on the \"field AGE of card 6\" area of forms viewer", () => clickArea(page, "field AGE of card 6", el("forms viewer")));
      await session.step(36, "Then the current column should be \"AGE\"", () => currentColumnIs(page, "AGE"));
      await session.step(37, "And \"USUBJID\" of the current row should be \"X0273T21001500015\"", () => currentRowValue(page, "USUBJID", "X0273T21001500015"));
      await session.step(38, "When user makes row 6 current", () => makeRowCurrent(page, 6));
      await session.step(39, "Then the \"current record\" reading of forms viewer should be 6", () => readingIs(page, "current record", el("forms viewer"), 6));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a header label sets the current column", async () => {
      await session.step(43, "When user clicks on the \"label HEIGHT\" area of forms viewer", () => clickArea(page, "label HEIGHT", el("forms viewer")));
      await session.step(44, "Then the current column should be \"HEIGHT\"", () => currentColumnIs(page, "HEIGHT"));
      await session.step(45, "And the \"current record\" reading of forms viewer should be 6", () => readingIs(page, "current record", el("forms viewer"), 6));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control-clicking a card toggles its row's selection", async () => {
      await session.step(49, "When user clicks on the \"current card\" area of forms viewer holding Control", () => clickAreaHolding(page, "current card", el("forms viewer"), "Control"));
      await session.step(50, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
      await session.step(51, "And no rows where \"USUBJID\" is \"X0273T21000500008\" should be selected", () => noneOfSelected(page, "USUBJID", "X0273T21000500008"));
      await session.step(52, "When user clicks on the \"current card\" area of forms viewer holding Control", () => clickAreaHolding(page, "current card", el("forms viewer"), "Control"));
      await session.step(53, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(54, "And all rows where \"USUBJID\" is \"X0273T21000500008\" should be selected", () => allOfSelected(page, "USUBJID", "X0273T21000500008"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift-clicking a card selects every row up to it", async () => {
      await session.step(58, "When user clicks on the \"current card\" area of forms viewer holding Shift", () => clickAreaHolding(page, "current card", el("forms viewer"), "Shift"));
      await session.step(59, "Then 6 rows should be selected", () => selectedRowCount(page, 6));
      await session.step(60, "And all rows where \"USUBJID\" is \"X0273T21000300003\" should be selected", () => allOfSelected(page, "USUBJID", "X0273T21000300003"));
      await session.step(61, "And no rows where \"USUBJID\" is \"X0273T21000900003\" should be selected", () => noneOfSelected(page, "USUBJID", "X0273T21000900003"));
      await session.step(62, "When user selects rows where \"USUBJID\" is one of \"X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015\"", () => selectWhereOneOf(page, "USUBJID", "X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015"));
      await session.step(63, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control-Shift-clicking a card clears every row up to it", async () => {
      await session.step(67, "When user clicks on the \"current card\" area of forms viewer holding Control and Shift", () => chordClick(page, "current card", "Control", "Shift"));
      await session.step(68, "Then 2 rows should be selected", () => selectedRowCount(page, 2));
      await session.step(69, "And all rows where \"USUBJID\" is \"X0273T21000900003\" should be selected", () => allOfSelected(page, "USUBJID", "X0273T21000900003"));
      await session.step(70, "And no rows where \"USUBJID\" is \"X0273T21000400001\" should be selected", () => noneOfSelected(page, "USUBJID", "X0273T21000400001"));
      await session.step(71, "When user selects rows where \"USUBJID\" is one of \"X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015\"", () => selectWhereOneOf(page, "USUBJID", "X0273T21000400001, X0273T21000500008, X0273T21000900003, X0273T21001500015"));
      await session.step(72, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering a card binds the mouse-over row and leaving releases it", async () => {
      await session.step(76, "When user hovers over the \"card 5\" area of forms viewer", () => hoverArea(page, "card 5", el("forms viewer")));
      await session.step(77, "Then the \"mouse-over record\" reading of forms viewer should be 12", () => readingIs(page, "mouse-over record", el("forms viewer"), 12));
      await session.step(78, "And the \"USUBJID of mouse-over card\" reading of forms viewer should be \"X0273T21000900003\"", () => readingReads(page, "USUBJID of mouse-over card", el("forms viewer"), "X0273T21000900003"));
      await session.step(79, "And the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000500008"));
      await session.step(80, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(81, "Then the \"mouse-over record\" reading of forms viewer should be \"\"", () => readingReads(page, "mouse-over record", el("forms viewer"), ""));
      await session.step(82, "And forms viewer should not have a \"field USUBJID of mouse-over card\" area", () => hasNoArea(page, el("forms viewer"), "field USUBJID of mouse-over card"));
      await session.step(83, "And forms viewer should have a \"current card\" area", () => hasArea(page, el("forms viewer"), "current card"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A hovered grid row fills the mouse-over card", async () => {
      await session.step(87, "When user hovers over the \"cell 12 of AGE\" area of grid", () => hoverArea(page, "cell 12 of AGE", el("grid")));
      await session.step(88, "Then the \"mouse-over record\" reading of forms viewer should be 12", () => readingIs(page, "mouse-over record", el("forms viewer"), 12));
      await session.step(89, "And the \"USUBJID of mouse-over card\" reading of forms viewer should be \"X0273T21000900003\"", () => readingReads(page, "USUBJID of mouse-over card", el("forms viewer"), "X0273T21000900003"));
      await session.step(90, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(91, "Then the \"mouse-over record\" reading of forms viewer should be \"\"", () => readingReads(page, "mouse-over record", el("forms viewer"), ""));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Mouse Over Row off drops the card and the grid hover with it", async () => {
      await session.step(95, "When user sets \"Show Mouse Over Row\" property of forms viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row", el("forms viewer"), "false"));
      await session.step(96, "Then the \"cards\" reading of forms viewer should be lower than before", () => readingLower(page, "cards", el("forms viewer")));
      await session.step(97, "And forms viewer should not have a \"mouse-over card\" area", () => hasNoArea(page, el("forms viewer"), "mouse-over card"));
      await session.step(98, "When user takes a snapshot of forms viewer", () => takeSnapshot(page, el("forms viewer")));
      await session.step(99, "And user hovers over the \"cell 12 of AGE\" area of grid", () => hoverArea(page, "cell 12 of AGE", el("grid")));
      await session.step(100, "Then the \"cards\" reading of forms viewer should be the same as before", () => readingSame(page, "cards", el("forms viewer")));
      await session.step(101, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(102, "And user sets \"Show Mouse Over Row\" property of forms viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row", el("forms viewer"), "true"));
      await session.step(103, "Then the \"cards\" reading of forms viewer should be higher than before", () => readingHigher(page, "cards", el("forms viewer")));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Current Row off drops the current card", async () => {
      await session.step(107, "When user sets \"Show Current Row\" property of forms viewer to \"false\"", () => setProperty(page, "Show Current Row", el("forms viewer"), "false"));
      await session.step(108, "Then forms viewer should not have a \"current card\" area", () => hasNoArea(page, el("forms viewer"), "current card"));
      await session.step(109, "And the \"cards\" reading of forms viewer should be lower than before", () => readingLower(page, "cards", el("forms viewer")));
      await session.step(110, "When user sets \"Show Current Row\" property of forms viewer to \"true\"", () => setProperty(page, "Show Current Row", el("forms viewer"), "true"));
      await session.step(111, "Then forms viewer should have a \"current card\" area", () => hasArea(page, el("forms viewer"), "current card"));
      await session.step(112, "And the \"cards\" reading of forms viewer should be higher than before", () => readingHigher(page, "cards", el("forms viewer")));
      await session.step(113, "And the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000500008\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000500008"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
