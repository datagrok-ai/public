/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/manage/collections.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.app.monomer-collections]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {collectionHolds, noCollectionOnServer, noSuchCollection} from '../../bindings/monomer-libs.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openApp} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Monomer collections", () => {
  const session = feature(test, "features/manage/collections.feature", import.meta.url);
  test("Monomer collections", {tag: ["@journey", "@realizes:bio.app.monomer-collections"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "And no \"bdd-test-collection\" monomer collection is on the server", () => noCollectionOnServer(page, "bdd-test-collection"));
    await session.step(12, "And user opens the \"Monomer Collections\" app", () => openApp(page, "Monomer Collections"));
    await session.step(13, "Then \"Canonical AAs\" card should be visible", () => shouldBe(page, el("\"Canonical AAs\" card"), "visible"));
    await session.step(14, "And \"New Collection\" card should be visible", () => shouldBe(page, el("\"New Collection\" card"), "visible"));
    await run.scenario("A new collection is made from the dialog and saved on the server", async () => {
      await session.step(17, "When user clicks on \"New Collection\" card", () => clickOn(page, el("\"New Collection\" card")));
      await session.step(18, "Then \"New Monomer Collection\" dialog should be visible", () => shouldBe(page, el("\"New Monomer Collection\" dialog"), "visible"));
      await session.step(19, "And \"Polymer Type\" input in \"New Monomer Collection\" dialog should have value \"PEPTIDE\"", () => shouldHaveValue(page, el("\"Polymer Type\" input in \"New Monomer Collection\" dialog"), "PEPTIDE"));
      await session.step(20, "When user types \"bdd-test-collection\" into Name input in \"New Monomer Collection\" dialog", () => typeInto(page, "bdd-test-collection", el("Name input in \"New Monomer Collection\" dialog")));
      await session.step(21, "And user types \"made by the bdd feature\" into Description input in \"New Monomer Collection\" dialog", () => typeInto(page, "made by the bdd feature", el("Description input in \"New Monomer Collection\" dialog")));
      await session.step(22, "And user types \"A\" into Search input in \"New Monomer Collection\" dialog", () => typeInto(page, "A", el("Search input in \"New Monomer Collection\" dialog")));
      await session.step(23, "And user presses Enter in Search input in \"New Monomer Collection\" dialog", () => pressKeyIn(page, "Enter", el("Search input in \"New Monomer Collection\" dialog")));
      await session.step(24, "And user types \"G\" into Search input in \"New Monomer Collection\" dialog", () => typeInto(page, "G", el("Search input in \"New Monomer Collection\" dialog")));
      await session.step(25, "And user presses Enter in Search input in \"New Monomer Collection\" dialog", () => pressKeyIn(page, "Enter", el("Search input in \"New Monomer Collection\" dialog")));
      await session.step(26, "Then \"New Monomer Collection\" dialog should contain text \"2 monomer(s) selected\"", () => shouldContainText(page, el("\"New Monomer Collection\" dialog"), "2 monomer(s) selected"));
      await session.step(27, "When user clicks on OK button in \"New Monomer Collection\" dialog", () => clickOn(page, el("OK button in \"New Monomer Collection\" dialog")));
      await session.step(28, "Then \"New Monomer Collection\" dialog should be hidden", () => shouldBe(page, el("\"New Monomer Collection\" dialog"), "hidden"));
      await session.step(29, "And \"bdd-test-collection\" card should become visible", () => shouldBe(page, el("\"bdd-test-collection\" card"), "visible"));
      await session.step(30, "And \"bdd-test-collection\" card should contain text \"2 monomer(s)\"", () => shouldContainText(page, el("\"bdd-test-collection\" card"), "2 monomer(s)"));
      await session.step(31, "And \"bdd-test-collection\" card should contain text \"made by the bdd feature\"", () => shouldContainText(page, el("\"bdd-test-collection\" card"), "made by the bdd feature"));
      await session.step(32, "And the \"bdd-test-collection\" monomer collection should hold monomers \"A, G\"", () => collectionHolds(page, "bdd-test-collection", "A, G"));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A card selects on click", async () => {
      await session.step(36, "When user clicks on \"bdd-test-collection\" card", () => clickOn(page, el("\"bdd-test-collection\" card")));
      await session.step(37, "Then \"bdd-test-collection\" card should be selected", () => shouldBe(page, el("\"bdd-test-collection\" card"), "selected"));
      await session.step(38, "And \"Canonical AAs\" card should not be selected", () => shouldNotBe(page, el("\"Canonical AAs\" card"), "selected"));
    });
    await run.scenario("Delete removes the collection after confirmation", async () => {
      await session.step(41, "When user clicks on Delete button in \"bdd-test-collection\" card", () => clickOn(page, el("Delete button in \"bdd-test-collection\" card")));
      await session.step(42, "Then \"Delete Collection\" dialog should be visible", () => shouldBe(page, el("\"Delete Collection\" dialog"), "visible"));
      await session.step(43, "And \"Delete Collection\" dialog should contain text \"bdd-test-collection\"", () => shouldContainText(page, el("\"Delete Collection\" dialog"), "bdd-test-collection"));
      await session.step(44, "When user clicks on OK button in \"Delete Collection\" dialog", () => clickOn(page, el("OK button in \"Delete Collection\" dialog")));
      await session.step(45, "Then \"bdd-test-collection\" card should become hidden", () => shouldBe(page, el("\"bdd-test-collection\" card"), "hidden"));
      await session.step(46, "And there should be no \"bdd-test-collection\" monomer collection on the server", () => noSuchCollection(page, "bdd-test-collection"));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
