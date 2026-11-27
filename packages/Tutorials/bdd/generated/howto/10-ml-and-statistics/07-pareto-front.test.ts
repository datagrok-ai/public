/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/10-ml-and-statistics/07-pareto-front.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hoverOver, selectIn, shouldBe, shouldHaveValue, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {maxInRow, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noBalloons, noErrors, painted, propertyShouldBe, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · ML", () => {
  const session = feature(test, "features/howto/10-ml-and-statistics/07-pareto-front.feature", import.meta.url);
  test("Find the best trade-offs with a Pareto front", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "When the video says \"Cheap or powerful: the Pareto front finds the best compromises\" at the bottom", () => videoSaysAt(page, "Cheap or powerful: the Pareto front finds the best compromises", "bottom"));
    await session.step(28, "And user picks \"ML > Pareto Front...\" from the top menu", () => pickFromTopMenu(page, "ML > Pareto Front..."), undefined, "Pareto Front is in the ML menu");
    await session.step(29, "Then pareto front viewer should be visible", () => shouldBe(page, el("pareto front viewer"), "visible"));
    await session.step(30, "Then \"price\" input should have value \"minimize\"", () => shouldHaveValue(page, el("\"price\" input"), "minimize"));
    await session.step(32, "When user selects \"\" in \"highway.mpg\" input", () => selectIn(page, "", el("\"highway.mpg\" input")), undefined, "Optimize: minimize or maximize each column; leave mileage out");
    await session.step(34, "And user selects \"maximize\" in horsepower input", () => selectIn(page, "maximize", el("horsepower input")), undefined, "Horsepower: higher is better. Price stays at minimize");
    await session.step(35, "Then \"Minimize\" property of pareto front viewer should be \"price\"", () => propertyShouldBe(page, "Minimize", el("pareto front viewer"), "price"));
    await session.step(36, "And \"Maximize\" property of pareto front viewer should be \"horsepower\"", () => propertyShouldBe(page, "Maximize", el("pareto front viewer"), "horsepower"));
    await session.step(37, "And scatter plot viewer in pareto front viewer should be painted", () => painted(page, el("scatter plot viewer in pareto front viewer")));
    await session.step(38, "When the video says \"Green: no other car is both cheaper and more powerful\" at the bottom", () => videoSaysAt(page, "Green: no other car is both cheaper and more powerful", "bottom"));
    await session.step(40, "And user hovers over the \"marker of row 12\" area of scatter plot viewer in pareto front viewer", () => hoverArea(page, "marker of row 12", el("scatter plot viewer in pareto front viewer")), undefined, "Grey Mercedes: 123 hp for 31,600");
    await session.step(41, "Then the \"hovered row\" reading of scatter plot viewer in pareto front viewer should be 12", () => readingIs(page, "hovered row", el("scatter plot viewer in pareto front viewer"), 12));
    await session.step(42, "And the value of \"model\" column in row 12 should be \"mercedes\"", () => valueInRow(page, "model", 12, "mercedes"));
    await session.step(43, "And the value of \"Pareto optimality\" column in row 12 should be \"non-optimal\"", () => valueInRow(page, "Pareto optimality", 12, "non-optimal"));
    await session.step(45, "When user hovers over the \"marker of row 1\" area of scatter plot viewer in pareto front viewer", () => hoverArea(page, "marker of row 1", el("scatter plot viewer in pareto front viewer")), undefined, "Alfa Romeo beats it: more power for about half the price");
    await session.step(46, "Then the \"hovered row\" reading of scatter plot viewer in pareto front viewer should be 1", () => readingIs(page, "hovered row", el("scatter plot viewer in pareto front viewer"), 1));
    await session.step(47, "And the value of \"model\" column in row 1 should be \"alfaromeo\"", () => valueInRow(page, "model", 1, "alfaromeo"));
    await session.step(48, "And the value of \"Pareto optimality\" column in row 1 should be \"optimal\"", () => valueInRow(page, "Pareto optimality", 1, "optimal"));
    await session.step(50, "When user hovers over the \"marker of row 19\" area of scatter plot viewer in pareto front viewer", () => hoverArea(page, "marker of row 19", el("scatter plot viewer in pareto front viewer")), undefined, "Porsche: the dearest, but no car has more power");
    await session.step(51, "Then the \"hovered row\" reading of scatter plot viewer in pareto front viewer should be 19", () => readingIs(page, "hovered row", el("scatter plot viewer in pareto front viewer"), 19));
    await session.step(52, "And the value of \"model\" column in row 19 should be \"porsche\"", () => valueInRow(page, "model", 19, "porsche"));
    await session.step(53, "And the value of \"Pareto optimality\" column in row 19 should be \"optimal\"", () => valueInRow(page, "Pareto optimality", 19, "optimal"));
    await session.step(54, "And \"horsepower\" column should have its maximum in row 19", () => maxInRow(page, "horsepower", 19));
    await session.step(55, "And \"price\" column should have its maximum in row 19", () => maxInRow(page, "price", 19));
    await session.step(57, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")), undefined, "Below, one line per car across its columns, green on the front");
    await session.step(58, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
