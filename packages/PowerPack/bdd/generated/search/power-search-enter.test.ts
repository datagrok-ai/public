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
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
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
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "When user types \"QA\" into home search", () => typeInto(page, "QA", el("home search")));
    await session.step(31, "Then the search suggestions should be \"PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "PDB ID, e.g. 4AKZ"));
    await session.step(32, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(33, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(34, "Then the search should have finished", () => searchFinished(page));
    await session.step(35, "And home search should have value \"QA\"", () => shouldHaveValue(page, el("home search"), "QA"));
    await session.step(36, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(37, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(38, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"new\" with its suggestions shown and none highlighted finds functions and help pages [query=new, suggestions=New Users Today | New users This Month | New users This Year | New users last 3 months | New user last 7 days | New users yesterday]", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "When user types \"new\" into home search", () => typeInto(page, "new", el("home search")));
    await session.step(31, "Then the search suggestions should be \"New Users Today | New users This Month | New users This Year | New users last 3 months | New user last 7 days | New users yesterday\"", () => suggestionsAre(page, "New Users Today | New users This Month | New users This Year | New users last 3 months | New user last 7 days | New users yesterday"));
    await session.step(32, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(33, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(34, "Then the search should have finished", () => searchFinished(page));
    await session.step(35, "And home search should have value \"new\"", () => shouldHaveValue(page, el("home search"), "new"));
    await session.step(36, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(37, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(38, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"user\" with its suggestions shown and none highlighted finds functions and help pages [query=user, suggestions=DGUSER-{User Login}]", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "When user types \"user\" into home search", () => typeInto(page, "user", el("home search")));
    await session.step(31, "Then the search suggestions should be \"DGUSER-{User Login}\"", () => suggestionsAre(page, "DGUSER-{User Login}"));
    await session.step(32, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(33, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(34, "Then the search should have finished", () => searchFinished(page));
    await session.step(35, "And home search should have value \"user\"", () => shouldHaveValue(page, el("home search"), "user"));
    await session.step(36, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(37, "And the page address should contain \"search?q=\"", () => urlShouldContain(page, "search?q="));
    await session.step(38, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"a\", with no suggestion shown, finds functions and help pages", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(49, "When user types \"a\" into home search", () => typeInto(page, "a", el("home search")));
    await session.step(50, "Then the search should have finished", () => searchFinished(page));
    await session.step(51, "And no search suggestion should be shown", () => noSuggestions(page));
    await session.step(52, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(53, "Then home search should have value \"a\"", () => shouldHaveValue(page, el("home search"), "a"));
    await session.step(54, "And the search results should list the categories \"Functions, Help\"", () => searchListsCategories(page, "Functions, Help"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"1+1\" with its suggestion shown and none highlighted finds nothing", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(59, "When user types \"1+1\" into home search", () => typeInto(page, "1+1", el("home search")));
    await session.step(60, "Then the search suggestions should be \"PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "PDB ID, e.g. 4AKZ"));
    await session.step(61, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(62, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(63, "Then the search should have finished", () => searchFinished(page));
    await session.step(64, "And home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(65, "And the search results should show nothing", () => searchShowsNothing(page));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
    await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on \"Project[0-9]+\", with no suggestion shown, finds nothing", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(70, "When user types \"Project[0-9]+\" into home search", () => typeInto(page, "Project[0-9]+", el("home search")));
    await session.step(71, "Then the search should have finished", () => searchFinished(page));
    await session.step(72, "And no search suggestion should be shown", () => noSuggestions(page));
    await session.step(73, "When user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(74, "Then home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
    await session.step(75, "And the search results should show nothing", () => searchShowsNothing(page));
    await session.step(76, "And no errors should have been logged", () => noErrors(page));
    await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("\"new\" lists the Add New Column function", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(80, "When user types \"new\" into home search", () => typeInto(page, "new", el("home search")));
    await session.step(81, "Then the search should have finished", () => searchFinished(page));
    await session.step(82, "And the \"Functions\" category of the search results should list \"Add New Column\"", () => categoryListsItem(page, "Functions", "Add New Column"));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
  test("The arrow keys walk the suggestions, and Enter takes the highlighted one", {tag: ["@realizes:powerpack.search.power-pack", "@realizes:powerpack.view.welcome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(86, "When user types \"DG\" into home search", () => typeInto(page, "DG", el("home search")));
    await session.step(87, "Then the search suggestions should be \"DGUSER-{User Login} | PDB ID, e.g. 4AKZ\"", () => suggestionsAre(page, "DGUSER-{User Login} | PDB ID, e.g. 4AKZ"));
    await session.step(88, "And the highlighted search suggestion should be \"none\"", () => suggestionHighlighted(page, "none"));
    await session.step(89, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(90, "Then the highlighted search suggestion should be \"DGUSER-{User Login}\"", () => suggestionHighlighted(page, "DGUSER-{User Login}"));
    await session.step(91, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(92, "Then the highlighted search suggestion should be \"PDB ID, e.g. 4AKZ\"", () => suggestionHighlighted(page, "PDB ID, e.g. 4AKZ"));
    await session.step(93, "When user presses ArrowUp in home search", () => pressKeyIn(page, "ArrowUp", el("home search")));
    await session.step(94, "Then the highlighted search suggestion should be \"DGUSER-{User Login}\"", () => suggestionHighlighted(page, "DGUSER-{User Login}"));
    await session.step(95, "When user presses ArrowDown in home search", () => pressKeyIn(page, "ArrowDown", el("home search")));
    await session.step(96, "And user presses Enter in home search", () => pressKeyIn(page, "Enter", el("home search")));
    await session.step(97, "Then home search should have value \"4AKZ\"", () => shouldHaveValue(page, el("home search"), "4AKZ"));
    await session.step(98, "And the page address should contain \"search?q=4AKZ\"", () => urlShouldContain(page, "search?q=4AKZ"));
    await session.step(99, "And the search should have finished", () => searchFinished(page));
    await session.step(100, "And \"PDB: 4AKZ\" link in home search results should be visible", () => shouldBe(page, el("\"PDB: 4AKZ\" link in home search results"), "visible"));
    await session.step(101, "And no errors should have been logged", () => noErrors(page));
    await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
