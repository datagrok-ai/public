/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/formula-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, GROK-17004, GROK-20931]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {acceptCompletion, highlightStandsOut, highlightsExactly, highlightsInColor} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, holdsCode, hoverText, pasteInto, pressKeyIn, shouldBe, shouldContainText, shouldNotBe, typeAtCaret, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoElement} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The formula editor of Add New Column on demog-1000: autocomplete, hints and column highlighting", () => {
  const session = feature(test, "features/add-new-column/formula-editor.feature", import.meta.url);
  test("The formula editor of Add New Column on demog-1000: autocomplete, hints and column highlighting", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@realizes:GROK-17004", "@realizes:GROK-20931"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(29, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
    await session.step(30, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await run.scenario("A typed letter offers the functions that start with it, and Enter inserts one", async () => {
      await session.step(33, "When user types \"a\" into formula editor", () => typeInto(page, "a", el("formula editor")));
      await session.step(34, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(35, "And \"Abs\" completion should be visible", () => shouldBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(36, "And \"Acos\" completion should be visible", () => shouldBe(page, el("\"Acos\" completion"), "visible"));
      await session.step(37, "And \"Avg\" completion should be visible", () => shouldBe(page, el("\"Avg\" completion"), "visible"));
      await session.step(38, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(39, "Then formula editor should hold the code \"Abs(x)\"", () => holdsCode(page, el("formula editor"), "Abs(x)"));
      await session.step(40, "And \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(41, "And completion list should be hidden", () => shouldBe(page, el("completion list"), "hidden"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on an offered function inserts it", async () => {
      await session.step(45, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(46, "And user types \"a\" into formula editor", () => typeInto(page, "a", el("formula editor")));
      await session.step(47, "Then \"Acos\" completion should be visible", () => shouldBe(page, el("\"Acos\" completion"), "visible"));
      await session.step(48, "When user clicks on \"Acos\" completion", () => clickOn(page, el("\"Acos\" completion")));
      await session.step(49, "Then formula editor should hold the code \"Acos(x)\"", () => holdsCode(page, el("formula editor"), "Acos(x)"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Space offers the functions on an empty editor", async () => {
      await session.step(53, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(54, "And user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(55, "Then completion list should be hidden", () => shouldBe(page, el("completion list"), "hidden"));
      await session.step(56, "When user presses ControlLeft+Space in formula editor", () => pressKeyIn(page, "ControlLeft+Space", el("formula editor")));
      await session.step(57, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(58, "And \"Abs\" completion should be visible", () => shouldBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"$\" offers the columns of the table, not functions", async () => {
      await session.step(62, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(63, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(64, "And user types \"$\" into formula editor", () => typeInto(page, "$", el("formula editor")));
      await session.step(65, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(66, "And \"HEIGHT\" completion should be visible", () => shouldBe(page, el("\"HEIGHT\" completion"), "visible"));
      await session.step(67, "And \"WEIGHT\" completion should be visible", () => shouldBe(page, el("\"WEIGHT\" completion"), "visible"));
      await session.step(68, "And \"AGE\" completion should be visible", () => shouldBe(page, el("\"AGE\" completion"), "visible"));
      await session.step(69, "And \"Abs\" completion should not be visible", () => shouldNotBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(70, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pointer over an inserted function shows its signature", async () => {
      await session.step(74, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(75, "And user types \"a\" into formula editor", () => typeInto(page, "a", el("formula editor")));
      await session.step(76, "Then \"Abs\" completion should be visible", () => shouldBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(77, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(78, "Then formula editor should hold the code \"Abs(x)\"", () => holdsCode(page, el("formula editor"), "Abs(x)"));
      await session.step(79, "And \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(80, "And formula hint should contain text \"Abs(x:\"", () => shouldContainText(page, el("formula hint"), "Abs(x:"));
      await session.step(81, "When user hovers over the text \"Abs\" in formula editor", () => hoverText(page, "Abs", el("formula editor")));
      await session.step(82, "Then signature tooltip should be visible", () => shouldBe(page, el("signature tooltip"), "visible"));
      await session.step(83, "And signature tooltip should contain text \"Abs(x:\"", () => shouldContainText(page, el("signature tooltip"), "Abs(x:"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pasted ${col} reference is highlighted, a bare name is not", async () => {
      await session.step(87, "When user pastes \"Abs(age)\" into formula editor", () => pasteInto(page, "Abs(age)", el("formula editor")));
      await session.step(88, "Then formula editor should hold the code \"Abs(age)\"", () => holdsCode(page, el("formula editor"), "Abs(age)"));
      await session.step(89, "And formula editor should highlight the column references \"\"", () => highlightsExactly(page, el("formula editor"), ""));
      await session.step(90, "When user pastes \"Abs(${age})\" into formula editor", () => pasteInto(page, "Abs(${age})", el("formula editor")));
      await session.step(91, "Then formula editor should hold the code \"Abs(${age})\"", () => holdsCode(page, el("formula editor"), "Abs(${age})"));
      await session.step(92, "And formula editor should highlight the column references \"${age}\"", () => highlightsExactly(page, el("formula editor"), "${age}"));
      await session.step(93, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(94, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pasted $[col] reference is highlighted", async () => {
      await session.step(98, "When user pastes \"Avg($[age])\" into formula editor", () => pasteInto(page, "Avg($[age])", el("formula editor")));
      await session.step(99, "Then formula editor should hold the code \"Avg($[age])\"", () => holdsCode(page, el("formula editor"), "Avg($[age])"));
      await session.step(100, "And formula editor should highlight the column references \"$[age]\"", () => highlightsExactly(page, el("formula editor"), "$[age]"));
      await session.step(101, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(102, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column picked from the \"${\" autocomplete is highlighted", async () => {
      await session.step(106, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(107, "And user types \"Round(\" into formula editor", () => typeInto(page, "Round(", el("formula editor")));
      await session.step(108, "And user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(109, "And user types \"${\" at the caret", () => typeAtCaret(page, "${"));
      await session.step(110, "Then \"HEIGHT\" completion should be visible", () => shouldBe(page, el("\"HEIGHT\" completion"), "visible"));
      await session.step(111, "When user clicks on \"HEIGHT\" completion", () => clickOn(page, el("\"HEIGHT\" completion")));
      await session.step(112, "Then formula editor should contain text \"Round(${HEIGHT}\"", () => shouldContainText(page, el("formula editor"), "Round(${HEIGHT}"));
      await session.step(113, "And formula editor should highlight the column references \"${HEIGHT}\"", () => highlightsExactly(page, el("formula editor"), "${HEIGHT}"));
      await session.step(114, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(115, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table parameter offers the open tables, and the column parameter after it that table's columns", async () => {
      await session.step(120, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(121, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(122, "And user types \"Lookup(\" into formula editor", () => typeInto(page, "Lookup(", el("formula editor")));
      await session.step(123, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(124, "And \"cars\" completion should be visible", () => shouldBe(page, el("\"cars\" completion"), "visible"));
      await session.step(125, "And \"demog-1000\" completion should be visible", () => shouldBe(page, el("\"demog-1000\" completion"), "visible"));
      await session.step(126, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(127, "Then formula editor should hold the code \"Lookup(\\\"cars\\\"\"", () => holdsCode(page, el("formula editor"), "Lookup(\"cars\""));
      await session.step(128, "When user types \", \" at the caret", () => typeAtCaret(page, ", "));
      await session.step(129, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(130, "And \"price\" completion should be visible", () => shouldBe(page, el("\"price\" completion"), "visible"));
      await session.step(131, "And \"model\" completion should be visible", () => shouldBe(page, el("\"model\" completion"), "visible"));
      await session.step(132, "And \"AGE\" completion should not be visible", () => shouldNotBe(page, el("\"AGE\" completion"), "visible"));
      await session.step(133, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(134, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column dropped from the grid's header is highlighted", async () => {
      await session.step(137, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(138, "And user types \"Sin(\" into formula editor", () => typeInto(page, "Sin(", el("formula editor")));
      await session.step(139, "And user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(140, "And user drags the \"header WEIGHT\" area of grid onto formula editor", () => dragAreaOntoElement(page, "header WEIGHT", el("grid"), el("formula editor")));
      await session.step(141, "Then formula editor should contain text \"Sin(${WEIGHT}\"", () => shouldContainText(page, el("formula editor"), "Sin(${WEIGHT}"));
      await session.step(142, "And formula editor should highlight the column references \"${WEIGHT}\"", () => highlightsExactly(page, el("formula editor"), "${WEIGHT}"));
      await session.step(143, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(144, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(145, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(146, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
