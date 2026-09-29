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
import {categoryListsItem, noSuggestions, searchFinished, searchListsCategories, searchShowsNothing, suggestionHighlighted, suggestionsAre} from '../../bindings/search.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKeyIn, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {urlShouldContain} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Enter in the Home search box", () => {
  const session = feature(test, "features/search/power-search-enter.feature", import.meta.url);
  test("Enter on \"QA\" with its suggestions shown and none highlighted finds functions and help pages [query=QA, suggestions=PDB ID, e.g. 4AKZ]", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "When user types \"QA\" into home search", () => typeInto(page, "QA", el("home search")));
    await session.step(33, "Then the search suggestions should be \"PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "PDB ID, e.g. 4AKZ"));
    await session.step(34, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(35, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(36, "Then the search should have finished", () => searchFinished(page));
    await session.step(37, "And home search should have value \"QA\"", () => shouldHaveValue(page, el("home search"), "QA"));
    await session.step(38, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(39, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(40, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"new\" with its suggestions shown and none highlighted finds functions and help pages [query=new, suggestions=New Users Today | New users This Month | New users This Year | New users last 3 months | New user last 7 days | New users yesterday]", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "When user types \"new\" into home search", () => typeInto(page, "new", el("home search")));
    await session.step(33, "Then the search suggestions should be \"New Users Today | New users This Month | New users This Year | New users last 3 months | New user last 7 days | New users yesterday\"", () => suggestionsAre(page, "New Users Today | New users This Month | New users This Year | New users last 3 months | New user last 7 days | New users yesterday"));
    await session.step(34, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(35, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(36, "Then the search should have finished", () => searchFinished(page));
    await session.step(37, "And home search should have value \"new\"", () => shouldHaveValue(page, el("home search"), "new"));
    await session.step(38, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(39, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(40, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"use\" with its suggestions shown and none highlighted finds functions and help pages [query=use, suggestions=DGUSER-{User Login}]", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "When user types \"use\" into home search", () => typeInto(page, "use", el("home search")));
    await session.step(33, "Then the search suggestions should be \"DGUSER-{User Login}\"", () => suggestionsAre(page, "DGUSER-{User Login}"));
    await session.step(34, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(35, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(36, "Then the search should have finished", () => searchFinished(page));
    await session.step(37, "And home search should have value \"use\"", () => shouldHaveValue(page, el("home search"), "use"));
    await session.step(38, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(39, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(40, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"a\" finds functions and help pages", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(51, "When user types \"a\" into home search", () => typeInto(page, "a", el("home search")));
    await session.step(52, "Then the search should have finished", () => searchFinished(page));
    await session.step(53, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(54, "Then home search should have value \"a\"", () => shouldHaveValue(page, el("home search"), "a"));
    await session.step(55, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
    await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"1+1\" with its suggestion shown and none highlighted finds nothing", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(60, "When user types \"1+1\" into home search", () => typeInto(page, "1+1", el("home search")));
    await session.step(61, "Then the search suggestions should be \"PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "PDB ID, e.g. 4AKZ"));
    await session.step(62, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(63, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(64, "Then the search should have finished", () => searchFinished(page));
    await session.step(65, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(66, "And the search results should show nothing", () => searchShowsNothing(page));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
    await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"Project[0-9]+\", with no suggestion shown, finishes without an error", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(71, "When user types \"Project[0-9]+\" into home search", () => typeInto(page, "Project[0-9]+", el("home search")));
    await session.step(72, "Then the search should have finished", () => searchFinished(page));
    await session.step(73, "And no search suggestion should be shown", () => noSuggestions(page));
    await session.step(74, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(75, "Then home search should have value \"Project[0-9]+\"", () => shouldHaveValue(page, el("home search"), "Project[0-9]+"));
    await session.step(76, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(77, "And the search should have finished", () => searchFinished(page));
    await session.step(78, "And no errors should have been logged", () => noErrors(page));
    await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("\"new\" lists the Add New Column function", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(82, "When user types \"new\" into home search", () => typeInto(page, "new", el("home search")));
    await session.step(83, "Then the search should have finished", () => searchFinished(page));
    await session.step(84, "And the \"Functions\" category of the search results should list \"Add New Column\"", () => categoryListsItem(page, "Functions", "Add New Column"));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
  test("The arrow keys walk the suggestions, and Enter takes the highlighted one", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(88, "When user types \"DG\" into home search", () => typeInto(page, "DG", el("home search")));
    await session.step(89, "Then the search suggestions should be \"DGUSER-{User Login} | PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "DGUSER-{User Login} | PDB ID, e.g. 4AKZ"));
    await session.step(90, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(91, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(92, "Then the highlighted search suggestion should be \"DGUSER-{User Login}\"", () => suggestionHighlighted(page, "DGUSER-{User Login}"));
    await session.step(93, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(94, "Then the highlighted search suggestion should be \"PDB ID, e.g. 4AKZ\"", () => suggestionHighlighted(page, "PDB ID, e.g. 4AKZ"));
    await session.step(95, "When user presses ArrowUp in home search", () => pressKeyIn(page, "ArrowUp", el("home search")));
    await session.step(96, "Then the highlighted search suggestion should be \"DGUSER-{User Login}\"", () => suggestionHighlighted(page, "DGUSER-{User Login}"));
    await session.step(97, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(98, "Then home search should have value \"DGUSER-\"", () => shouldHaveValue(page, el("home search"), "DGUSER-"));
    await session.step(99, "And the page address should contain \"search?q=DGUSER-\"", () => urlShouldContain(page, "search?q=DGUSER-"));
    await session.step(100, "And the search should have finished", () => searchFinished(page));
    await session.step(101, "And no errors should have been logged", () => noErrors(page));
    await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
