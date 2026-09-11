/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/pareto-front-objectives.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [eda.viewer.pareto-front]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnIncomplete, columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCalculated} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, propertyShouldBe, propertyShouldContain, propertyShouldNotBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pareto front objectives", () => {
  const session = feature(test, "features/pareto-front-objectives.feature", import.meta.url);
  test("Pareto front objectives", {tag: ["@journey", "@eda", "@realizes:eda.viewer.pareto-front", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user picks \"ML > Pareto Front...\" from the top menu", () => pickFromTopMenu(page, "ML > Pareto Front..."));
    await session.step(29, "Then pareto front viewer should be visible", () => shouldBe(page, el("pareto front viewer"), "visible"));
    await run.scenario("The viewer's properties come in their categories, and only numeric columns are offered", async () => {
      await session.step(32, "When user clicks on settings icon of pareto front viewer", () => clickOn(page, el("settings icon of pareto front viewer")));
      await session.step(33, "Then \"Objectives\" category in context panel should be visible", () => shouldBe(page, el("\"Objectives\" category in context panel"), "visible"));
      await session.step(34, "And \"Axes\" category in context panel should be visible", () => shouldBe(page, el("\"Axes\" category in context panel"), "visible"));
      await session.step(35, "And \"Labels\" category in context panel should be visible", () => shouldBe(page, el("\"Labels\" category in context panel"), "visible"));
      await session.step(36, "And \"Legend\" category in context panel should be visible", () => shouldBe(page, el("\"Legend\" category in context panel"), "visible"));
      await session.step(37, "And \"Description\" category in context panel should be visible", () => shouldBe(page, el("\"Description\" category in context panel"), "visible"));
      await session.step(38, "Given \"Objectives\" category in context panel is expanded", () => isExpanded(page, el("\"Objectives\" category in context panel")));
      await session.step(39, "Then \"Minimize\" property in context panel should contain text \"2 / 16\"", () => shouldContainText(page, el("\"Minimize\" property in context panel"), "2 / 16"));
      await session.step(40, "And \"Maximize\" property in context panel should contain text \"0 / 16\"", () => shouldContainText(page, el("\"Maximize\" property in context panel"), "0 / 16"));
      await session.step(41, "When user clicks on \"...\" button in \"Maximize\" property in context panel", () => clickOn(page, el("\"...\" button in \"Maximize\" property in context panel")));
      await session.step(42, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(43, "And the \"rows\" reading of grid viewer in \"Select columns...\" dialog should be 16", () => readingIs(page, "rows", el("grid viewer in \"Select columns...\" dialog"), 16));
      await session.step(44, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"diesel\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "diesel"));
      await session.step(45, "When user clicks on CANCEL button in \"Select columns...\" dialog", () => clickOn(page, el("CANCEL button in \"Select columns...\" dialog")));
      await session.step(46, "Then \"Maximize\" property of pareto front viewer should be \"\"", () => propertyShouldBe(page, "Maximize", el("pareto front viewer"), ""));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Maximizing what is minimized is refused with a warning, and choosing again clears it", async () => {
      await session.step(50, "When user clicks on \"...\" button in \"Maximize\" property in context panel", () => clickOn(page, el("\"...\" button in \"Maximize\" property in context panel")));
      await session.step(51, "And user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(52, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(53, "Then \"Maximize\" property in context panel should contain text \"16 / 16\"", () => shouldContainText(page, el("\"Maximize\" property in context panel"), "16 / 16"));
      await session.step(54, "And \"Maximize\" property of pareto front viewer should contain \"price\"", () => propertyShouldContain(page, "Maximize", el("pareto front viewer"), "price"));
      await session.step(55, "And pareto front viewer should contain text \"Cannot minimize and maximize features at the same time\"", () => shouldContainText(page, el("pareto front viewer"), "Cannot minimize and maximize features at the same time"));
      await session.step(56, "And pareto front viewer should contain text \"highway.mpg\"", () => shouldContainText(page, el("pareto front viewer"), "highway.mpg"));
      await session.step(57, "When user clicks on \"...\" button in \"Maximize\" property in context panel", () => clickOn(page, el("\"...\" button in \"Maximize\" property in context panel")));
      await session.step(58, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(59, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(60, "Then \"Maximize\" property of pareto front viewer should be \"\"", () => propertyShouldBe(page, "Maximize", el("pareto front viewer"), ""));
      await session.step(61, "And pareto front viewer should not contain text \"Cannot minimize and maximize\"", () => shouldNotContainText(page, el("pareto front viewer"), "Cannot minimize and maximize"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An axis and the labels chosen by hand turn their automatic choice off", async () => {
      await session.step(65, "Given \"Axes\" category in context panel is expanded", () => isExpanded(page, el("\"Axes\" category in context panel")));
      await session.step(66, "When user selects \"horsepower\" in \"X Axis\" property in context panel", () => selectIn(page, "horsepower", el("\"X Axis\" property in context panel")));
      await session.step(67, "Then \"X Axis\" property of pareto front viewer should be \"horsepower\"", () => propertyShouldBe(page, "X Axis", el("pareto front viewer"), "horsepower"));
      await session.step(68, "And \"Auto Axes Selection\" property of pareto front viewer should not be \"true\"", () => propertyShouldNotBe(page, "Auto Axes Selection", el("pareto front viewer"), "true"));
      await session.step(69, "Given \"Labels\" category in context panel is expanded", () => isExpanded(page, el("\"Labels\" category in context panel")));
      await session.step(70, "When user clicks on \"...\" button in \"Label Columns\" property in context panel", () => clickOn(page, el("\"...\" button in \"Label Columns\" property in context panel")));
      await session.step(71, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(72, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(73, "Then \"Label Columns\" property of pareto front viewer should be \"\"", () => propertyShouldBe(page, "Label Columns", el("pareto front viewer"), ""));
      await session.step(74, "And \"Auto Labels Selection\" property of pareto front viewer should not be \"true\"", () => propertyShouldNotBe(page, "Auto Labels Selection", el("pareto front viewer"), "true"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An empty column is not offered as an objective", async () => {
      await session.step(79, "When user adds a calculated column \"empty\" with formula \"If(true, null, 0)\"", () => addCalculated(page, "empty", "If(true, null, 0)"));
      await session.step(80, "Then \"empty\" column should have type \"int\"", () => columnType(page, "empty", "int"));
      await session.step(81, "And \"empty\" column should have missing values", () => columnIncomplete(page, "empty"));
      await session.step(82, "When user clicks on \"...\" button in \"Maximize\" property in context panel", () => clickOn(page, el("\"...\" button in \"Maximize\" property in context panel")));
      await session.step(83, "Then the \"rows\" reading of grid viewer in \"Select columns...\" dialog should be 16", () => readingIs(page, "rows", el("grid viewer in \"Select columns...\" dialog"), 16));
    }, {knownFailure: true});
    run.finish();
  });
});
