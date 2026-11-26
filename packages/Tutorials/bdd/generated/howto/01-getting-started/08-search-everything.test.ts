/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/08-search-everything.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
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
import {homeSearchCategories, homeSearchFinished, homeSearchTableLoaded} from '../../../bindings/howto-basics.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, hoverOver, shouldBe, typeInto, typeKeyByKey, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {clickViewTab, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/08-search-everything.feature", import.meta.url);
  test("Search everything from Home", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And simple mode is off", () => simpleModeOff(page));
    await session.step(23, "When the video says \"The Home search box looks through the whole server\"", () => videoSays(page, "The Home search box looks through the whole server"));
    await session.step(24, "And user clicks on home search", () => clickOn(page, el("home search")));
    await session.step(26, "And user types \"scatter plot\" into home search", () => typeInto(page, "scatter plot", el("home search")), undefined, "Type what you look for: results come by kind");
    await session.step(27, "Then the Home search should have finished", () => homeSearchFinished(page));
    await session.step(28, "And the Home search results should list the kinds \"Functions, Demos, Help, Samples\"", () => homeSearchCategories(page, "Functions, Demos, Help, Samples"));
    await session.step(30, "When user hovers over \"Scatterplot\" link inside home search results", () => hoverOver(page, el("\"Scatterplot\" link inside home search results")), undefined, "Functions, demos, code samples and help pages");
    await session.step(32, "And user types \"Threshold\" key by key into home search", () => typeKeyByKey(page, "Threshold", el("home search")), undefined, "A word from a name finds dashboards, files and queries");
    await session.step(33, "Then the Home search should have finished", () => homeSearchFinished(page));
    await session.step(34, "And the Home search results should list the kinds \"Projects, Files, Queries\"", () => homeSearchCategories(page, "Projects, Files, Queries"));
    await session.step(36, "When user hovers over second \"ADCs by IC50 Threshold\" link inside home search results", () => hoverOver(page, el("second \"ADCs by IC50 Threshold\" link inside home search results")), undefined, "Point at a dashboard to preview it");
    await session.step(38, "And user double-clicks on second \"ADCs by IC50 Threshold\" link inside home search results", () => doubleClickOn(page, el("second \"ADCs by IC50 Threshold\" link inside home search results")), undefined, "Double-click to open it");
    await session.step(39, "Then the \"ADCs by IC50 Threshold\" view should be current", () => viewIsCurrent(page, "ADCs by IC50 Threshold"));
    await session.step(41, "When user clicks on the tab of the \"Home\" view", () => clickViewTab(page, "Home"), undefined, "The house tab takes you back to Home");
    await session.step(42, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(44, "When user types \"biologics\" key by key into home search", () => typeKeyByKey(page, "biologics", el("home search")), undefined, "Some queries answer a sentence. Start typing one");
    await session.step(46, "And user hovers over \"biologics assays for ${organism}\" menu item inside context menu", () => hoverOver(page, el("\"biologics assays for ${organism}\" menu item inside context menu")), undefined, "The box suggests the sentences it knows");
    await session.step(48, "And user types \"biologics assays for E. coli\" key by key into home search", () => typeKeyByKey(page, "biologics assays for E. coli", el("home search")), undefined, "Finish the sentence and the query runs");
    await session.step(49, "Then the Home search should have finished", () => homeSearchFinished(page));
    await session.step(50, "And the table the Home search answered with should have loaded", () => homeSearchTableLoaded(page));
    await session.step(52, "When user hovers over home search results", () => hoverOver(page, el("home search results")), undefined, "Its answer: the assays run in E. coli");
    await session.step(54, "And user types \"GROKADC-000001\" into home search", () => typeInto(page, "GROKADC-000001", el("home search")), undefined, "Type an ID to see its record");
    await session.step(55, "Then the Home search should have finished", () => homeSearchFinished(page));
    await session.step(56, "When user hovers over home search results", () => hoverOver(page, el("home search results")));
    await session.step(58, "And user types \"GROKMOL-000001, GROKMOL-000002, GROKMOL-000003\" into home search", () => typeInto(page, "GROKMOL-000001, GROKMOL-000002, GROKMOL-000003", el("home search")), undefined, "Several IDs, separated by commas, show a card each");
    await session.step(59, "Then the Home search should have finished", () => homeSearchFinished(page));
    await session.step(60, "When user hovers over home search results", () => hoverOver(page, el("home search results")));
    await session.step(62, "And user clicks on \"Open table\" link inside home search results", () => clickOn(page, el("\"Open table\" link inside home search results")), undefined, "Open table puts them in a table");
    await session.step(63, "Then the table should have 3 rows", () => rowCount(page, 3));
    await session.step(65, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The IDs stay links: click one to see its record");
    await session.step(66, "And user clicks on the \"cell 2 of DG_BIOLOGICS_DRUG_ID\" area of grid", () => clickArea(page, "cell 2 of DG_BIOLOGICS_DRUG_ID", el("grid")));
    await session.step(67, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
});
