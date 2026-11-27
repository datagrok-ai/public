/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/08-end-to-end/02-explore-a-new-dataset.feature
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
import {clickOn, expand, hoverText, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, onlyOfAnySelected, removeColumn, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelOpen, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, dragSelectionBetweenAreas, hoverArea, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · End to end", () => {
  const session = feature(test, "features/howto/08-end-to-end/02-explore-a-new-dataset.feature", import.meta.url);
  test("Explore a new dataset", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(21, "And user removes \"diesel\" column", () => removeColumn(page, "diesel"));
    await session.step(22, "And user removes \"two.doors\" column", () => removeColumn(page, "two.doors"));
    await session.step(23, "And user removes \"hatchback\" column", () => removeColumn(page, "hatchback"));
    await session.step(24, "And user removes \"wheel.base\" column", () => removeColumn(page, "wheel.base"));
    await session.step(25, "And user removes \"length\" column", () => removeColumn(page, "length"));
    await session.step(26, "And user removes \"width\" column", () => removeColumn(page, "width"));
    await session.step(27, "And user removes \"height\" column", () => removeColumn(page, "height"));
    await session.step(28, "And user removes \"peak.rpm\" column", () => removeColumn(page, "peak.rpm"));
    await session.step(29, "And user removes \"symbol\" column", () => removeColumn(page, "symbol"));
    await session.step(30, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(31, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(32, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(33, "When the video says \"A table you have not seen before. Get to know it\"", () => videoSays(page, "A table you have not seen before. Get to know it"));
    await session.step(35, "And user hovers over the text \"Columns: 8\" in status bar", () => hoverText(page, "Columns: 8", el("status bar")), undefined, "The status bar: 30 rows, 8 columns");
    await session.step(37, "And user clicks on the \"header price\" area of grid", () => clickArea(page, "header price", el("grid")), undefined, "Click a column: Stats gives its range and average");
    await session.step(38, "Then the context panel should show \"price\"", () => contextPanelShows(page, "price"));
    await session.step(39, "When user expands Stats pane in context panel", () => expand(page, el("Stats pane in context panel")));
    await session.step(40, "Then Stats pane in context panel should contain text \"avg\"", () => shouldContainText(page, el("Stats pane in context panel"), "avg"));
    await session.step(41, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")));
    await session.step(43, "And user clicks on \"histogram\" icon in toolbox", () => clickOn(page, el("\"histogram\" icon in toolbox")), undefined, "A histogram shows the spread: most cars are cheap");
    await session.step(44, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(45, "When user picks \"price\" in the \"value\" column selector of histogram viewer", () => pickInColumnSelector(page, "price", "value", el("histogram viewer")));
    await session.step(47, "And user drags a selection box from the \"bin 17\" area to the \"bin 20\" area of histogram viewer", () => dragSelectionBetweenAreas(page, "bin 17", "bin 20", el("histogram viewer")), undefined, "Three cost far more. Shift+drag over their bars to select them");
    await session.step(48, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
    await session.step(49, "And only rows where \"model\" is one of \"jaguar, mercedes, porsche\" should be selected", () => onlyOfAnySelected(page, "model", "jaguar, mercedes, porsche"));
    await session.step(50, "When the video says \"The grid marks the selected rows too\"", () => videoSays(page, "The grid marks the selected rows too"));
    await session.step(51, "And user hovers over the \"cell 10 of model\" area of grid", () => hoverArea(page, "cell 10 of model", el("grid")));
    await session.step(53, "And user clicks on correlation plot icon on toolbox", () => clickOn(page, el("correlation plot icon on toolbox")), undefined, "Why so expensive: a correlation plot shows what goes with price");
    await session.step(54, "Then correlation plot viewer should be visible", () => shouldBe(page, el("correlation plot viewer"), "visible"));
    await session.step(56, "When user hovers over the \"cell city.mpg x price\" area of correlation plot viewer", () => hoverArea(page, "cell city.mpg x price", el("correlation plot viewer")), undefined, "Blue: one falls as the other rises, as mileage does with price");
    await session.step(58, "And user hovers over the \"cell eng.size x price\" area of correlation plot viewer", () => hoverArea(page, "cell eng.size x price", el("correlation plot viewer")), undefined, "Red: they rise together. Engine size goes with price most");
    await session.step(60, "And user double-clicks on the \"cell eng.size x price\" area of correlation plot viewer", () => doubleClickArea(page, "cell eng.size x price", el("correlation plot viewer")), undefined, "Double-click the cell to plot the two");
    await session.step(61, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(62, "When the video says \"The three selected cars have the biggest engines: top right\"", () => videoSays(page, "The three selected cars have the biggest engines: top right"));
    await session.step(63, "And user hovers over the \"marker of row 19\" area of scatter plot viewer", () => hoverArea(page, "marker of row 19", el("scatter plot viewer")));
    await session.step(65, "And user picks \"Select > Selection to Filter\" from the top menu", () => pickFromTopMenu(page, "Select > Selection to Filter"), undefined, "Selection to Filter shows only those rows, everywhere");
    await session.step(66, "Then 3 rows should pass the filter", () => filterPasses(page, 3));
    await session.step(67, "And every selected row should pass the filter", () => selectedPassFilter(page));
    await session.step(68, "When the video says \"Jaguar, Mercedes and Porsche: the priciest, with the biggest engines\"", () => videoSays(page, "Jaguar, Mercedes and Porsche: the priciest, with the biggest engines"));
    await session.step(69, "And user hovers over the \"cell 10 of model\" area of grid", () => hoverArea(page, "cell 10 of model", el("grid")));
    await session.step(70, "Then no errors should have been logged", () => noErrors(page));
  });
});
