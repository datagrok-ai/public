/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/formula-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, GROK-17004]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {acceptCompletion, dragAreaOnto, highlightStandsOut, highlightsExactly, highlightsInColor, holdsFormula, hoverText, typeAtCaret} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, pasteInto, pressKeyIn, shouldBe, shouldContainText, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The formula editor of Add New Column on demog: autocomplete, hints and column highlighting", () => {
  const session = feature(test, "features/add-new-column/formula-editor.feature", import.meta.url);
  test("The formula editor of Add New Column on demog: autocomplete, hints and column highlighting", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@realizes:GROK-17004"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
    await session.step(20, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await run.scenario("A typed letter offers the functions that start with it, and Enter inserts one", async () => {
      await session.step(23, "When user types \"a\" into formula editor", () => typeInto(page, "a", el("formula editor")));
      await session.step(24, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(25, "And \"Abs\" completion should be visible", () => shouldBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(26, "And \"Acos\" completion should be visible", () => shouldBe(page, el("\"Acos\" completion"), "visible"));
      await session.step(27, "And \"Avg\" completion should be visible", () => shouldBe(page, el("\"Avg\" completion"), "visible"));
      await session.step(28, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(29, "Then formula editor should hold the formula \"Abs(x)\"", () => holdsFormula(page, el("formula editor"), "Abs(x)"));
      await session.step(30, "And \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(31, "And completion list should be hidden", () => shouldBe(page, el("completion list"), "hidden"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on an offered function inserts it", async () => {
      await session.step(35, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(36, "And user types \"a\" into formula editor", () => typeInto(page, "a", el("formula editor")));
      await session.step(37, "Then \"Acos\" completion should be visible", () => shouldBe(page, el("\"Acos\" completion"), "visible"));
      await session.step(38, "When user clicks on \"Acos\" completion", () => clickOn(page, el("\"Acos\" completion")));
      await session.step(39, "Then formula editor should hold the formula \"Acos(x)\"", () => holdsFormula(page, el("formula editor"), "Acos(x)"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Space offers the functions on an empty editor", async () => {
      await session.step(43, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(44, "And user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(45, "Then completion list should be hidden", () => shouldBe(page, el("completion list"), "hidden"));
      await session.step(46, "When user presses Control+Space in formula editor", () => pressKeyIn(page, "Control+Space", el("formula editor")));
      await session.step(47, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(48, "And \"Abs\" completion should be visible", () => shouldBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"$\" offers the columns of the table, not functions", async () => {
      await session.step(52, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(53, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(54, "And user types \"$\" into formula editor", () => typeInto(page, "$", el("formula editor")));
      await session.step(55, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(56, "And \"HEIGHT\" completion should be visible", () => shouldBe(page, el("\"HEIGHT\" completion"), "visible"));
      await session.step(57, "And \"WEIGHT\" completion should be visible", () => shouldBe(page, el("\"WEIGHT\" completion"), "visible"));
      await session.step(58, "And \"AGE\" completion should be visible", () => shouldBe(page, el("\"AGE\" completion"), "visible"));
      await session.step(59, "And \"Abs\" completion should not be visible", () => shouldNotBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(60, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pointer over an inserted function shows its signature", async () => {
      await session.step(64, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(65, "And user types \"a\" into formula editor", () => typeInto(page, "a", el("formula editor")));
      await session.step(66, "Then \"Abs\" completion should be visible", () => shouldBe(page, el("\"Abs\" completion"), "visible"));
      await session.step(67, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(68, "Then formula editor should hold the formula \"Abs(x)\"", () => holdsFormula(page, el("formula editor"), "Abs(x)"));
      await session.step(69, "And \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(70, "And formula hint should contain text \"Abs(x:\"", () => shouldContainText(page, el("formula hint"), "Abs(x:"));
      await session.step(71, "When user hovers over the text \"Abs\" in formula editor", () => hoverText(page, "Abs", el("formula editor")));
      await session.step(72, "Then signature tooltip should be visible", () => shouldBe(page, el("signature tooltip"), "visible"));
      await session.step(73, "And signature tooltip should contain text \"Abs(x:\"", () => shouldContainText(page, el("signature tooltip"), "Abs(x:"));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pasted ${col} reference is highlighted, a bare name is not", async () => {
      await session.step(77, "When user pastes \"Abs(age)\" into formula editor", () => pasteInto(page, "Abs(age)", el("formula editor")));
      await session.step(78, "Then formula editor should hold the formula \"Abs(age)\"", () => holdsFormula(page, el("formula editor"), "Abs(age)"));
      await session.step(79, "And formula editor should highlight the column references \"\"", () => highlightsExactly(page, el("formula editor"), ""));
      await session.step(80, "When user pastes \"Abs(${age})\" into formula editor", () => pasteInto(page, "Abs(${age})", el("formula editor")));
      await session.step(81, "Then formula editor should hold the formula \"Abs(${age})\"", () => holdsFormula(page, el("formula editor"), "Abs(${age})"));
      await session.step(82, "And formula editor should highlight the column references \"${age}\"", () => highlightsExactly(page, el("formula editor"), "${age}"));
      await session.step(83, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(84, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pasted $[col] reference is highlighted", async () => {
      await session.step(88, "When user pastes \"Avg($[age])\" into formula editor", () => pasteInto(page, "Avg($[age])", el("formula editor")));
      await session.step(89, "Then formula editor should hold the formula \"Avg($[age])\"", () => holdsFormula(page, el("formula editor"), "Avg($[age])"));
      await session.step(90, "And formula editor should highlight the column references \"$[age]\"", () => highlightsExactly(page, el("formula editor"), "$[age]"));
      await session.step(91, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(92, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column picked from the \"${\" autocomplete is highlighted", async () => {
      await session.step(96, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(97, "And user types \"Round(\" into formula editor", () => typeInto(page, "Round(", el("formula editor")));
      await session.step(98, "And user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(99, "And user types \"${\" at the caret", () => typeAtCaret(page, "${"));
      await session.step(100, "Then \"HEIGHT\" completion should be visible", () => shouldBe(page, el("\"HEIGHT\" completion"), "visible"));
      await session.step(101, "When user clicks on \"HEIGHT\" completion", () => clickOn(page, el("\"HEIGHT\" completion")));
      await session.step(102, "Then formula editor should contain text \"Round(${HEIGHT}\"", () => shouldContainText(page, el("formula editor"), "Round(${HEIGHT}"));
      await session.step(103, "And formula editor should highlight the column references \"${HEIGHT}\"", () => highlightsExactly(page, el("formula editor"), "${HEIGHT}"));
      await session.step(104, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(105, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column dropped from the grid's header is highlighted", async () => {
      await session.step(109, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(110, "And user types \"Sin(\" into formula editor", () => typeInto(page, "Sin(", el("formula editor")));
      await session.step(111, "And user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(112, "And user drags the \"header WEIGHT\" area of grid onto formula editor", () => dragAreaOnto(page, "header WEIGHT", el("grid"), el("formula editor")));
      await session.step(113, "Then formula editor should contain text \"Sin(${WEIGHT}\"", () => shouldContainText(page, el("formula editor"), "Sin(${WEIGHT}"));
      await session.step(114, "And formula editor should highlight the column references \"${WEIGHT}\"", () => highlightsExactly(page, el("formula editor"), "${WEIGHT}"));
      await session.step(115, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(116, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(117, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(118, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
