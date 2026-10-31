/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/edit-column-formula.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, public#434]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {previewComputes, previewLacks} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, holdsCode, isExpanded, pressKey, pressKeyIn, replaceCode, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, columnIsCurrentObject} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, doubleClickArea, noBalloons, noErrors, readingAsRemembered, readingNotAsRemembered, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Edit Column Formula: a calculated column renamed and given a new formula in the dialog, and recalculated by a grid edit", () => {
  const session = feature(test, "features/add-new-column/edit-column-formula.feature", import.meta.url);
  test("Edit Column Formula: a calculated column renamed and given a new formula in the dialog, and recalculated by a grid edit", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@realizes:public#434"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens a table \"calc\" with:", () => openTableOf(page, "calc", [["x","s"],["1","a"],["2","b"],["4","c"]]), [["x","s"],["1","a"],["2","b"],["4","c"]]);
    await session.step(30, "And user adds a calculated column \"dbl\" with formula \"${x} * 2\"", () => addCalculated(page, "dbl", "${x} * 2"));
    await session.step(31, "And user adds a calculated column \"next\" with formula \"${dbl} + 1\"", () => addCalculated(page, "next", "${dbl} + 1"));
    await session.step(32, "And user adds a calculated column \"run\" with formula \"CumSum(${x})\"", () => addCalculated(page, "run", "CumSum(${x})"));
    await session.step(33, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The Formula pane opens the dialog with the column's name and formula", async () => {
      await session.step(36, "When user clicks on the \"header dbl\" area of grid", () => clickArea(page, "header dbl", el("grid")));
      await session.step(37, "Then the context panel should show \"dbl\"", () => contextPanelShows(page, "dbl"));
      await session.step(38, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(39, "Then formula pane editor in Formula pane in context panel should hold the code \"${x} * 2\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${x} * 2"));
      await session.step(40, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(41, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(42, "And column name input should have value \"dbl\"", () => shouldHaveValue(page, el("column name input"), "dbl"));
      await session.step(43, "And formula editor should hold the code \"${x} * 2\"", () => holdsCode(page, el("formula editor"), "${x} * 2"));
      await session.step(44, "When user clicks on CANCEL button in \"Edit Column Formula\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Column Formula\" dialog")));
      await session.step(45, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog of the column over it starts from that column's formula, and CANCEL drops an edit", async () => {
      await session.step(49, "Given the \"next\" column is the current object", () => columnIsCurrentObject(page, "next"));
      await session.step(50, "Then the context panel should show \"next\"", () => contextPanelShows(page, "next"));
      await session.step(51, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(52, "Then formula pane editor in Formula pane in context panel should hold the code \"${dbl} + 1\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${dbl} + 1"));
      await session.step(53, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(54, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(55, "And column name input should have value \"next\"", () => shouldHaveValue(page, el("column name input"), "next"));
      await session.step(56, "And formula editor should hold the code \"${dbl} + 1\"", () => holdsCode(page, el("formula editor"), "${dbl} + 1"));
      await session.step(57, "When user replaces the code of formula editor with \"${x}\"", () => replaceCode(page, el("formula editor"), "${x}"));
      await session.step(58, "And user clicks on CANCEL button in \"Edit Column Formula\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Column Formula\" dialog")));
      await session.step(59, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(60, "And \"next\" column should have tag \"formula\" equal to \"${dbl} + 1\"", () => columnTag(page, "next", "formula", "${dbl} + 1"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A new name and a new formula rename and recalculate the column, and the column over it follows", async () => {
      await session.step(64, "Given the \"dbl\" column is the current object", () => columnIsCurrentObject(page, "dbl"));
      await session.step(65, "Then the context panel should show \"dbl\"", () => contextPanelShows(page, "dbl"));
      await session.step(66, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(67, "Then formula pane editor in Formula pane in context panel should hold the code \"${x} * 2\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${x} * 2"));
      await session.step(68, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(69, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(70, "When user types \"triple\" into column name input", () => typeInto(page, "triple", el("column name input")));
      await session.step(71, "And user replaces the code of formula editor with \"${x} * 3\"", () => replaceCode(page, el("formula editor"), "${x} * 3"));
      await session.step(72, "And user clicks on OK button in \"Edit Column Formula\" dialog", () => clickOn(page, el("OK button in \"Edit Column Formula\" dialog")));
      await session.step(73, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(74, "And the table should have a column \"triple\"", () => hasColumn(page, "triple"));
      await session.step(75, "And the table should not have a column \"dbl\"", () => hasNoColumn(page, "dbl"));
      await session.step(76, "And \"triple\" column should have tag \"formula\" equal to \"${x} * 3\"", () => columnTag(page, "triple", "formula", "${x} * 3"));
      await session.step(77, "And the value of \"triple\" column in row 3 should be \"12\"", () => valueInRow(page, "triple", 3, "12"));
      await session.step(78, "And \"next\" column should have tag \"formula\" equal to \"${triple} + 1\"", () => columnTag(page, "next", "formula", "${triple} + 1"));
      await session.step(79, "And the value of \"next\" column in row 3 should be \"13\"", () => valueInRow(page, "next", 3, "13"));
      await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A value edited in the grid recalculates the running total", async () => {
      await session.step(84, "Then the value of \"run\" column in row 3 should be \"7\"", () => valueInRow(page, "run", 3, "7"));
      await session.step(85, "When user double-clicks on the \"cell 1 of x\" area of grid", () => doubleClickArea(page, "cell 1 of x", el("grid")));
      await session.step(86, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(87, "And user types \"10\" into cell editor", () => typeInto(page, "10", el("cell editor")));
      await session.step(88, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(89, "Then the value of \"x\" column in row 1 should be \"10\"", () => valueInRow(page, "x", 1, "10"));
      await session.step(90, "And the value of \"run\" column in row 1 should be \"10\"", () => valueInRow(page, "run", 1, "10"));
      await session.step(91, "And the value of \"run\" column in row 3 should be \"16\"", () => valueInRow(page, "run", 3, "16"));
      await session.step(92, "And the value of \"triple\" column in row 1 should be \"30\"", () => valueInRow(page, "triple", 1, "30"));
      await session.step(93, "And the value of \"next\" column in row 1 should be \"31\"", () => valueInRow(page, "next", 1, "31"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A scatter plot on a calculated column redraws after its formula is changed in the dialog, and keeps the new values (public#434)", async () => {
      await session.step(98, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","x"],["Y","triple"]]), [["X","x"],["Y","triple"]]);
      await session.step(101, "When user remembers the \"y axis max\" reading of scatter plot viewer", () => rememberReading(page, "y axis max", el("scatter plot viewer")));
      await session.step(102, "Given the \"triple\" column is the current object", () => columnIsCurrentObject(page, "triple"));
      await session.step(103, "Then the context panel should show \"triple\"", () => contextPanelShows(page, "triple"));
      await session.step(104, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(105, "Then formula pane editor in Formula pane in context panel should hold the code \"${x} * 3\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${x} * 3"));
      await session.step(106, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(107, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(108, "When user replaces the code of formula editor with \"${x} * 30\"", () => replaceCode(page, el("formula editor"), "${x} * 30"));
      await session.step(109, "And user clicks on OK button in \"Edit Column Formula\" dialog", () => clickOn(page, el("OK button in \"Edit Column Formula\" dialog")));
      await session.step(110, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(111, "And the \"y axis max\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "y axis max", el("scatter plot viewer")));
      await session.step(112, "When user remembers the \"y axis max\" reading of scatter plot viewer", () => rememberReading(page, "y axis max", el("scatter plot viewer")));
      await session.step(113, "And user sets \"Y\" property of scatter plot viewer to \"x\"", () => setProperty(page, "Y", el("scatter plot viewer"), "x"));
      await session.step(114, "Then the \"y axis max\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "y axis max", el("scatter plot viewer")));
      await session.step(115, "When user sets \"Y\" property of scatter plot viewer to \"triple\"", () => setProperty(page, "Y", el("scatter plot viewer"), "triple"));
      await session.step(116, "Then the \"y axis max\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "y axis max", el("scatter plot viewer")));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The preview of a column made in the dialog shows it under its own name", async () => {
      await session.step(120, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(121, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(122, "When user types \"half\" into column name input", () => typeInto(page, "half", el("column name input")));
      await session.step(123, "And user types \"${x} / 2\" into formula editor", () => typeInto(page, "${x} / 2", el("formula editor")));
      await session.step(124, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(125, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(126, "And the table should have a column \"half\"", () => hasColumn(page, "half"));
      await session.step(127, "Given the \"half\" column is the current object", () => columnIsCurrentObject(page, "half"));
      await session.step(128, "Then the context panel should show \"half\"", () => contextPanelShows(page, "half"));
      await session.step(129, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(130, "Then formula pane editor in Formula pane in context panel should hold the code \"${x} / 2\"", () => holdsCode(page, el("formula pane editor in Formula pane in context panel"), "${x} / 2"));
      await session.step(131, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
      await session.step(132, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
      await session.step(133, "And the preview grid should not show \"s\"", () => previewLacks(page, "s"));
      await session.step(134, "And the preview grid should not show \"half (2)\"", () => previewLacks(page, "half (2)"));
      await session.step(135, "And the preview grid should show \"half\" computed as numbers", () => previewComputes(page, "half"));
      await session.step(136, "When user clicks on CANCEL button in \"Edit Column Formula\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Column Formula\" dialog")));
      await session.step(137, "Then \"Edit Column Formula\" dialog should be hidden", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "hidden"));
      await session.step(138, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
