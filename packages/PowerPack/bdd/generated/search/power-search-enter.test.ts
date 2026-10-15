/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/power-search-enter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.search.power-pack, powerpack.view.welcome]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noSuggestions, searchFinished, searchListsCategories, suggestionHighlighted, suggestionsAre} from '../../bindings/search.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKeyIn, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {urlShouldContain} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Enter in the Home search box", () => {
  const session = feature(test, "features/search/power-search-enter.feature", import.meta.url);
  test("Enter on \"QA\" with its suggestions shown and none highlighted finds functions and help pages [query=QA, suggestions=PDB ID, e.g. 4AKZ]", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "When user types \"QA\" into home search", () => typeInto(page, "QA", el("home search")));
    await session.step(34, "Then the search suggestions should be \"PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "PDB ID, e.g. 4AKZ"));
    await session.step(35, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(36, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(37, "Then the search should have finished", () => searchFinished(page));
    await session.step(38, "And home search should have value \"QA\"", () => shouldHaveValue(page, el("home search"), "QA"));
    await session.step(39, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(40, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(41, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"a\" finds functions and help pages", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(50, "When user types \"a\" into home search", () => typeInto(page, "a", el("home search")));
    await session.step(51, "Then the search should have finished", () => searchFinished(page));
    await session.step(52, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(53, "Then home search should have value \"a\"", () => shouldHaveValue(page, el("home search"), "a"));
    await session.step(54, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"Project[0-9]+\", with no suggestion shown, finishes without an error", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(59, "When user types \"Project[0-9]+\" into home search", () => typeInto(page, "Project[0-9]+", el("home search")));
    await session.step(60, "Then the search should have finished", () => searchFinished(page));
    await session.step(61, "And no search suggestion should be shown", () => noSuggestions(page));
    await session.step(62, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(63, "Then home search should have value \"Project[0-9]+\"", () => shouldHaveValue(page, el("home search"), "Project[0-9]+"));
    await session.step(64, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(65, "And the search should have finished", () => searchFinished(page));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
    await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The arrow keys walk the suggestions, and Enter takes the highlighted one", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(70, "When user types \"DG\" into home search", () => typeInto(page, "DG", el("home search")));
    await session.step(71, "Then the search suggestions should be \"DGUSER-{User Login} | PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "DGUSER-{User Login} | PDB ID, e.g. 4AKZ"));
    await session.step(72, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(73, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(74, "Then the highlighted search suggestion should be \"DGUSER-{User Login}\"", () => suggestionHighlighted(page, "DGUSER-{User Login}"));
    await session.step(75, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(76, "Then the highlighted search suggestion should be \"PDB ID, e.g. 4AKZ\"", () => suggestionHighlighted(page, "PDB ID, e.g. 4AKZ"));
    await session.step(77, "When user presses ArrowUp in home search", () => pressKeyIn(page, "ArrowUp", el("home search")));
    await session.step(78, "Then the highlighted search suggestion should be \"DGUSER-{User Login}\"", () => suggestionHighlighted(page, "DGUSER-{User Login}"));
    await session.step(79, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(80, "Then home search should have value \"DGUSER-\"", () => shouldHaveValue(page, el("home search"), "DGUSER-"));
    await session.step(81, "And the page address should contain \"search?q=DGUSER-\"", () => urlShouldContain(page, "search?q=DGUSER-"));
    await session.step(82, "And the search should have finished", () => searchFinished(page));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
    await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
