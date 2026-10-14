/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/edit-column-formula.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, holdsCode, isExpanded, pressKey, pressKeyIn, replaceCode, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, columnIsCurrentObject} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Edit Column Formula: a calculated column renamed and given a new formula in the dialog, and recalculated by a grid edit", () => {
  const session = feature(test, "features/add-new-column/edit-column-formula.feature", import.meta.url);
  test("Edit Column Formula: a calculated column renamed and given a new formula in the dialog, and recalculated by a grid edit", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens a table \"calc\" with:", () => openTableOf(page, "calc", [["x","s"],["1","a"],["2","b"],["4","c"]]), [["x","s"],["1","a"],["2","b"],["4","c"]]);
    await session.step(27, "And user adds a calculated column \"dbl\" with formula \"${x} * 2\"", () => addCalculated(page, "dbl", "${x} * 2"));
    await session.step(28, "And user adds a calculated column \"next\" with formula \"${dbl} + 1\"", () => addCalculated(page, "next", "${dbl} + 1"));
    await session.step(29, "And user adds a calculated column \"run\" with formula \"CumSum(${x})\"", () => addCalculated(page, "run", "CumSum(${x})"));
    await session.step(30, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The Formula pane opens the dialog with the column's name and formula", async () => {
      await session.step(33, "When user clicks on the \"header dbl\" area of grid", () => clickArea(page, "header dbl", el("grid")));
      await session.step(34, "Then the context panel should show \"dbl\"", () => contextPanelShows(page, "dbl"));
      await session.step(35, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(36, "Then formula pane editor in Formula pane in context panel should hold the code \"${x} * 2\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${x} * 2"));
      await session.step(37, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(38, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(39, "And column name input should have value \"dbl\"", () => shouldHaveValue(page, el("column name input"), "dbl"));
      await session.step(40, "And formula editor should hold the code \"${x} * 2\"", () => holdsCode(page, el("formula editor"), "${x} * 2"));
      await session.step(41, "When user clicks on CANCEL button in \"Edit Column Formula\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Column Formula\" dialog")));
      await session.step(42, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog of the column over it starts from that column's formula, and CANCEL drops an edit", async () => {
      await session.step(46, "Given the \"next\" column is the current object", () => columnIsCurrentObject(page, "next"));
      await session.step(47, "Then the context panel should show \"next\"", () => contextPanelShows(page, "next"));
      await session.step(48, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(49, "Then formula pane editor in Formula pane in context panel should hold the code \"${dbl} + 1\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${dbl} + 1"));
      await session.step(50, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(51, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(52, "And column name input should have value \"next\"", () => shouldHaveValue(page, el("column name input"), "next"));
      await session.step(53, "And formula editor should hold the code \"${dbl} + 1\"", () => holdsCode(page, el("formula editor"), "${dbl} + 1"));
      await session.step(54, "When user replaces the code of formula editor with \"${x}\"", () => replaceCode(page, el("formula editor"), "${x}"));
      await session.step(55, "And user clicks on CANCEL button in \"Edit Column Formula\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Column Formula\" dialog")));
      await session.step(56, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(57, "And \"next\" column should have tag \"formula\" equal to \"${dbl} + 1\"", () => columnTag(page, "next", "formula", "${dbl} + 1"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A new name and a new formula rename and recalculate the column, and the column over it follows", async () => {
      await session.step(61, "Given the \"dbl\" column is the current object", () => columnIsCurrentObject(page, "dbl"));
      await session.step(62, "Then the context panel should show \"dbl\"", () => contextPanelShows(page, "dbl"));
      await session.step(63, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(64, "Then formula pane editor in Formula pane in context panel should hold the code \"${x} * 2\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${x} * 2"));
      await session.step(65, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(66, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(67, "When user types \"triple\" into column name input", () => typeInto(page, "triple", el("column name input")));
      await session.step(68, "And user replaces the code of formula editor with \"${x} * 3\"", () => replaceCode(page, el("formula editor"), "${x} * 3"));
      await session.step(69, "And user clicks on OK button in \"Edit Column Formula\" dialog", () => clickOn(page, el("OK button in \"Edit Column Formula\" dialog")));
      await session.step(70, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(71, "And the table should have a column \"triple\"", () => hasColumn(page, "triple"));
      await session.step(72, "And the table should not have a column \"dbl\"", () => hasNoColumn(page, "dbl"));
      await session.step(73, "And \"triple\" column should have tag \"formula\" equal to \"${x} * 3\"", () => columnTag(page, "triple", "formula", "${x} * 3"));
      await session.step(74, "And the value of \"triple\" column in row 3 should be \"12\"", () => valueInRow(page, "triple", 3, "12"));
      await session.step(75, "And \"next\" column should have tag \"formula\" equal to \"${triple} + 1\"", () => columnTag(page, "next", "formula", "${triple} + 1"));
      await session.step(76, "And the value of \"next\" column in row 3 should be \"13\"", () => valueInRow(page, "next", 3, "13"));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A value edited in the grid recalculates the running total", async () => {
      await session.step(81, "Then the value of \"run\" column in row 3 should be \"7\"", () => valueInRow(page, "run", 3, "7"));
      await session.step(82, "When user double-clicks on the \"cell 1 of x\" area of grid", () => doubleClickArea(page, "cell 1 of x", el("grid")));
      await session.step(83, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(84, "And user types \"10\" into cell editor", () => typeInto(page, "10", el("cell editor")));
      await session.step(85, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(86, "Then the value of \"x\" column in row 1 should be \"10\"", () => valueInRow(page, "x", 1, "10"));
      await session.step(87, "And the value of \"run\" column in row 1 should be \"10\"", () => valueInRow(page, "run", 1, "10"));
      await session.step(88, "And the value of \"run\" column in row 3 should be \"16\"", () => valueInRow(page, "run", 3, "16"));
      await session.step(89, "And the value of \"triple\" column in row 1 should be \"30\"", () => valueInRow(page, "triple", 1, "30"));
      await session.step(90, "And the value of \"next\" column in row 1 should be \"31\"", () => valueInRow(page, "next", 1, "31"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
