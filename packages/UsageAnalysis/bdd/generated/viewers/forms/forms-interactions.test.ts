/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-interactions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {recordCardRows} from '../../../bindings/forms.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, currentColumnIs, rowsRangeSelected, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, hoverArea, noErrors, pointerAway, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer mouse interactions and row binding", () => {
  const session = feature(test, "features/viewers/forms/forms-interactions.feature", import.meta.url);
  test("Forms viewer mouse interactions and row binding", {tag: ["@journey", "@viewers", "@realizes:viewers.forms", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(19, "And user makes row 13 current", () => makeRowCurrent(page, 13));
    await session.step(20, "And user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
    await session.step(21, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(22, "And the \"cards\" reading of forms viewer should be 7", () => readingIs(page, "cards", el("forms viewer"), 7));
    await session.step(23, "And the record cards of forms viewer should show rows \"215, 304, 428, 430, 512\"", () => recordCardRows(page, el("forms viewer"), "215, 304, 428, 430, 512"));
    await session.step(24, "And the \"record of card 1\" reading of forms viewer should be 13", () => readingIs(page, "record of card 1", el("forms viewer"), 13));
    await run.scenario("Clicking a card makes its row current", async () => {
      await session.step(27, "Then row 13 should be current", () => currentRowIs(page, 13));
      await session.step(28, "When user clicks on the \"card 4\" area of forms viewer", () => clickArea(page, "card 4", el("forms viewer")));
      await session.step(29, "Then row 304 should be current", () => currentRowIs(page, 304));
      await session.step(30, "And the \"record of card 1\" reading of forms viewer should be 304", () => readingIs(page, "record of card 1", el("forms viewer"), 304));
      await session.step(31, "And the \"SEX of card 1\" reading of forms viewer should be \"M\"", () => readingReads(page, "SEX of card 1", el("forms viewer"), "M"));
      await session.step(32, "And the \"AGE of card 1\" reading of forms viewer should be \"59\"", () => readingReads(page, "AGE of card 1", el("forms viewer"), "59"));
      await session.step(33, "When user makes row 13 current", () => makeRowCurrent(page, 13));
      await session.step(34, "Then the \"record of card 1\" reading of forms viewer should be 13", () => readingIs(page, "record of card 1", el("forms viewer"), 13));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a field makes that column and that row the current cell", async () => {
      await session.step(38, "When user clicks on the \"field AGE of card 5\" area of forms viewer", () => clickArea(page, "field AGE of card 5", el("forms viewer")));
      await session.step(39, "Then row 428 should be current", () => currentRowIs(page, 428));
      await session.step(40, "And the current column should be \"AGE\"", () => currentColumnIs(page, "AGE"));
      await session.step(41, "When user clicks on the \"field WEIGHT of card 6\" area of forms viewer", () => clickArea(page, "field WEIGHT of card 6", el("forms viewer")));
      await session.step(42, "Then row 430 should be current", () => currentRowIs(page, 430));
      await session.step(43, "And the current column should be \"WEIGHT\"", () => currentColumnIs(page, "WEIGHT"));
      await session.step(44, "When user makes row 13 current", () => makeRowCurrent(page, 13));
      await session.step(45, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a header label makes that column current without moving the row", async () => {
      await session.step(48, "When user clicks on the \"label HEIGHT\" area of forms viewer", () => clickArea(page, "label HEIGHT", el("forms viewer")));
      await session.step(49, "Then the current column should be \"HEIGHT\"", () => currentColumnIs(page, "HEIGHT"));
      await session.step(50, "And row 13 should be current", () => currentRowIs(page, 13));
      await session.step(51, "When user clicks on the \"label SEX\" area of forms viewer", () => clickArea(page, "label SEX", el("forms viewer")));
      await session.step(52, "Then the current column should be \"SEX\"", () => currentColumnIs(page, "SEX"));
      await session.step(53, "And row 13 should be current", () => currentRowIs(page, 13));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl and a click toggle that card's row out of the selection and back in", async () => {
      await session.step(57, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(58, "When user clicks on the \"card 3\" area of forms viewer holding Control", () => clickAreaHolding(page, "card 3", el("forms viewer"), "Control"));
      await session.step(59, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(60, "And the record cards of forms viewer should show rows \"304, 428, 430, 512\"", () => recordCardRows(page, el("forms viewer"), "304, 428, 430, 512"));
      await session.step(61, "When user clicks on the \"card 1\" area of forms viewer holding Control", () => clickAreaHolding(page, "card 1", el("forms viewer"), "Control"));
      await session.step(62, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift and a click select the run of rows up to the clicked card", async () => {
      await session.step(66, "Given user clears the row selection", () => clearSelection(page));
      await session.step(67, "And user makes row 13 current", () => makeRowCurrent(page, 13));
      await session.step(68, "When user clicks on the \"card 1\" area of forms viewer holding Shift", () => clickAreaHolding(page, "card 1", el("forms viewer"), "Shift"));
      await session.step(69, "Then 13 rows should be selected", () => selectedRowCount(page, 13));
      await session.step(70, "And rows 1 to 13 should be selected", () => rowsRangeSelected(page, 1, 13));
      await session.step(71, "When user clears the row selection", () => clearSelection(page));
      await session.step(72, "And user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
      await session.step(73, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl and Shift and a click clear the run instead of selecting it", async () => {
      await session.step(77, "Given user clears the row selection", () => clearSelection(page));
      await session.step(78, "And user makes row 13 current", () => makeRowCurrent(page, 13));
      await session.step(79, "When user clicks on the \"card 1\" area of forms viewer holding Shift", () => clickAreaHolding(page, "card 1", el("forms viewer"), "Shift"));
      await session.step(80, "Then rows 1 to 13 should be selected", () => rowsRangeSelected(page, 1, 13));
      await session.step(81, "When user makes row 9 current", () => makeRowCurrent(page, 9));
      await session.step(82, "And user clicks on the \"card 1\" area of forms viewer holding Control+Shift", () => clickAreaHolding(page, "card 1", el("forms viewer"), "Control+Shift"));
      await session.step(83, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(84, "And rows 10 to 13 should be selected", () => rowsRangeSelected(page, 10, 13));
      await session.step(85, "When user clears the row selection", () => clearSelection(page));
      await session.step(86, "And user makes row 13 current", () => makeRowCurrent(page, 13));
      await session.step(87, "And user selects rows where \"SEVERITY\" is \"Critical\"", () => selectWhereIs(page, "SEVERITY", "Critical"));
      await session.step(88, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(89, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering a card fills the mouse-over card, and leaving empties it again", async () => {
      await session.step(92, "Then the \"mouse-over record\" reading of forms viewer should be \"\"", () => readingReads(page, "mouse-over record", el("forms viewer"), ""));
      await session.step(93, "And the \"record of card 2\" reading of forms viewer should be \"\"", () => readingReads(page, "record of card 2", el("forms viewer"), ""));
      await session.step(94, "When user hovers over the \"card 4\" area of forms viewer", () => hoverArea(page, "card 4", el("forms viewer")));
      await session.step(95, "Then the \"mouse-over record\" reading of forms viewer should be 304", () => readingIs(page, "mouse-over record", el("forms viewer"), 304));
      await session.step(96, "And the \"record of card 2\" reading of forms viewer should be 304", () => readingIs(page, "record of card 2", el("forms viewer"), 304));
      await session.step(97, "And the \"USUBJID of mouse-over card\" reading of forms viewer should be \"X0273T29012500105\"", () => readingReads(page, "USUBJID of mouse-over card", el("forms viewer"), "X0273T29012500105"));
      await session.step(98, "And the \"card kind of card 2\" reading of forms viewer should be \"mouse-over\"", () => readingReads(page, "card kind of card 2", el("forms viewer"), "mouse-over"));
      await session.step(99, "When user hovers over the \"card 6\" area of forms viewer", () => hoverArea(page, "card 6", el("forms viewer")));
      await session.step(100, "Then the \"mouse-over record\" reading of forms viewer should be 430", () => readingIs(page, "mouse-over record", el("forms viewer"), 430));
      await session.step(101, "And the \"USUBJID of mouse-over card\" reading of forms viewer should be \"X0273T37001500013\"", () => readingReads(page, "USUBJID of mouse-over card", el("forms viewer"), "X0273T37001500013"));
      await session.step(102, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(103, "Then the \"mouse-over record\" reading of forms viewer should be \"\"", () => readingReads(page, "mouse-over record", el("forms viewer"), ""));
      await session.step(104, "And the \"record of card 2\" reading of forms viewer should be \"\"", () => readingReads(page, "record of card 2", el("forms viewer"), ""));
      await session.step(105, "And the \"cards\" reading of forms viewer should be 7", () => readingIs(page, "cards", el("forms viewer"), 7));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Mouse Over Row off takes the blank second card away", async () => {
      await session.step(109, "Then the \"cards\" reading of forms viewer should be 7", () => readingIs(page, "cards", el("forms viewer"), 7));
      await session.step(110, "And the \"card kind of card 2\" reading of forms viewer should be \"mouse-over\"", () => readingReads(page, "card kind of card 2", el("forms viewer"), "mouse-over"));
      await session.step(111, "When user sets \"showMouseOverRow\" property of forms viewer to \"false\"", () => setProperty(page, "showMouseOverRow", el("forms viewer"), "false"));
      await session.step(112, "Then the \"cards\" reading of forms viewer should be 6", () => readingIs(page, "cards", el("forms viewer"), 6));
      await session.step(113, "And the \"card kind of card 2\" reading of forms viewer should be \"record\"", () => readingReads(page, "card kind of card 2", el("forms viewer"), "record"));
      await session.step(114, "And the \"record of card 2\" reading of forms viewer should be 215", () => readingIs(page, "record of card 2", el("forms viewer"), 215));
      await session.step(115, "And the record cards of forms viewer should show rows \"215, 304, 428, 430, 512\"", () => recordCardRows(page, el("forms viewer"), "215, 304, 428, 430, 512"));
      await session.step(116, "When user hovers over the \"card 3\" area of forms viewer", () => hoverArea(page, "card 3", el("forms viewer")));
      await session.step(117, "Then the \"cards\" reading of forms viewer should be 6", () => readingIs(page, "cards", el("forms viewer"), 6));
      await session.step(118, "And the \"record of card 1\" reading of forms viewer should be 13", () => readingIs(page, "record of card 1", el("forms viewer"), 13));
      await session.step(119, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(120, "And user sets \"showMouseOverRow\" property of forms viewer to \"true\"", () => setProperty(page, "showMouseOverRow", el("forms viewer"), "true"));
      await session.step(121, "Then the \"cards\" reading of forms viewer should be 7", () => readingIs(page, "cards", el("forms viewer"), 7));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Current Row off leaves the record cards where they were", async () => {
      await session.step(132, "Given user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(133, "Then the \"mouse-over record\" reading of forms viewer should be \"\"", () => readingReads(page, "mouse-over record", el("forms viewer"), ""));
      await session.step(134, "And the \"cards\" reading of forms viewer should be 7", () => readingIs(page, "cards", el("forms viewer"), 7));
      await session.step(135, "When user sets \"showCurrentRow\" property of forms viewer to \"false\"", () => setProperty(page, "showCurrentRow", el("forms viewer"), "false"));
      await session.step(136, "Then the \"cards\" reading of forms viewer should be 6", () => readingIs(page, "cards", el("forms viewer"), 6));
      await session.step(137, "And the record cards of forms viewer should show rows \"215, 304, 428, 430, 512\"", () => recordCardRows(page, el("forms viewer"), "215, 304, 428, 430, 512"));
    }, {knownFailure: true});
    run.finish();
  });
});
