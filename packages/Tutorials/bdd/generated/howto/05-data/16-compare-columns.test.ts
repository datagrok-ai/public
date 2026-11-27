/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/16-compare-columns.feature
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
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {removeColumn, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/16-compare-columns.feature", import.meta.url);
  test("Compare two columns", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens tennis-players dataset", () => openDataset(page, ds("tennis-players")));
    await session.step(19, "And user removes \"image\" column", () => removeColumn(page, "image"));
    await session.step(20, "And user removes \"abstract\" column", () => removeColumn(page, "abstract"));
    await session.step(21, "And the grid shows only the columns \"name, country, birthPlace\"", () => gridShowsColumns(page, "name, country, birthPlace"));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "When user hovers over the \"cell 2 of birthPlace\" area of grid", () => hoverArea(page, "cell 2 of birthPlace", el("grid")), undefined, "country: who they play for; birthPlace: where they were born");
    await session.step(28, "And user picks \"Data > Compare Columns...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Columns..."), undefined, "Compare Columns: where two columns agree, row by row");
    await session.step(29, "Then Compare Columns dialog should be visible", () => shouldBe(page, el("Compare Columns dialog"), "visible"));
    await session.step(30, "When user selects \"tennis_players\" in \"Tables\" input in Compare Columns dialog", () => selectIn(page, "tennis_players", el("\"Tables\" input in Compare Columns dialog")));
    await session.step(32, "And user selects \"country\" in \"Columns\" input in Compare Columns dialog", () => selectIn(page, "country", el("\"Columns\" input in Compare Columns dialog")), undefined, "The two columns to compare");
    await session.step(33, "And user selects \"birthPlace\" in 4th input in Compare Columns dialog", () => selectIn(page, "birthPlace", el("4th input in Compare Columns dialog")));
    await session.step(34, "Then Compare Columns dialog should contain text \"Matched: 685\"", () => shouldContainText(page, el("Compare Columns dialog"), "Matched: 685"));
    await session.step(35, "And Compare Columns dialog should contain text \"Mismatched: 129\"", () => shouldContainText(page, el("Compare Columns dialog"), "Mismatched: 129"));
    await session.step(37, "When user clicks on \"Select\" link in Compare Columns dialog", () => clickOn(page, el("\"Select\" link in Compare Columns dialog")), undefined, "Select: the rows where they differ");
    await session.step(38, "Then 129 rows should be selected", () => selectedRowCount(page, 129));
    await session.step(39, "When user clicks on \"Close\" icon in Compare Columns dialog", () => clickOn(page, el("\"Close\" icon in Compare Columns dialog")));
    await session.step(40, "Then Compare Columns dialog should be hidden", () => shouldBe(page, el("Compare Columns dialog"), "hidden"));
    await session.step(41, "And the value of \"birthPlace\" column in row 2 should be \"Czech Republic\"", () => valueInRow(page, "birthPlace", 2, "Czech Republic"));
    await session.step(42, "And the value of \"country\" column in row 5 should be \"United Kingdom,United Kingdom\"", () => valueInRow(page, "country", 5, "United Kingdom,United Kingdom"));
    await session.step(44, "When user hovers over the \"cell 2 of birthPlace\" area of grid", () => hoverArea(page, "cell 2 of birthPlace", el("grid")), undefined, "Some players were born in another country");
    await session.step(46, "And user hovers over the \"cell 5 of country\" area of grid", () => hoverArea(page, "cell 5 of country", el("grid")), undefined, "Others are messy data: a country written twice");
    await session.step(47, "Then no errors should have been logged", () => noErrors(page));
  });
});
