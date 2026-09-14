/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, doubleClickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, noSpaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Searching spaces", () => {
  const session = feature(test, "features/spaces/spaces-search.feature", import.meta.url);
  test("Searching spaces", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(13, "And no space named \"BDD-Find, BDD-Miss, BDD-Find-Child\" is on the server", () => noSpaceOnServer(page, "BDD-Find, BDD-Miss, BDD-Find-Child"));
    await run.scenario("Two spaces to search among", async () => {
      await session.step(16, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(17, "And user enters \"BDD-Find\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Find", el("Name input in Create Space dialog")));
      await session.step(18, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(19, "Then 1 space named \"BDD-Find\" should be on the server", () => spacesOnServer(page, 1, "BDD-Find"));
      await session.step(20, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(21, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(22, "And user enters \"BDD-Miss\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Miss", el("Name input in Create Space dialog")));
      await session.step(23, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(24, "Then 1 space named \"BDD-Miss\" should be on the server", () => spacesOnServer(page, 1, "BDD-Miss"));
      await session.step(25, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    });
    await run.scenario("The Spaces list shows both", async () => {
      await session.step(28, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
      await session.step(29, "And user clicks on \"Refresh\" icon", () => clickOn(page, el("\"Refresh\" icon")));
      await session.step(30, "Then BDD-Find link in space gallery should be visible", () => shouldBe(page, el("BDD-Find link in space gallery"), "visible"));
      await session.step(31, "And BDD-Miss link in space gallery should be visible", () => shouldBe(page, el("BDD-Miss link in space gallery"), "visible"));
    });
    await run.scenario("A whole name keeps only that space", async () => {
      await session.step(34, "When user enters \"BDD-Find\" into space search", () => enterInto(page, "BDD-Find", el("space search")));
      await session.step(35, "Then BDD-Find link in space gallery should be visible", () => shouldBe(page, el("BDD-Find link in space gallery"), "visible"));
      await session.step(36, "And BDD-Miss link in space gallery should be absent", () => shouldBe(page, el("BDD-Miss link in space gallery"), "absent"));
    });
    await run.scenario("Part of a name still matches", async () => {
      await session.step(39, "When user enters \"BDD-Fi\" into space search", () => enterInto(page, "BDD-Fi", el("space search")));
      await session.step(40, "Then BDD-Find link in space gallery should be visible", () => shouldBe(page, el("BDD-Find link in space gallery"), "visible"));
      await session.step(41, "And BDD-Miss link in space gallery should be absent", () => shouldBe(page, el("BDD-Miss link in space gallery"), "absent"));
    });
    await run.scenario("A name nothing carries empties the list", async () => {
      await session.step(44, "When user enters \"zzz-no-such-space\" into space search", () => enterInto(page, "zzz-no-such-space", el("space search")));
      await session.step(45, "Then BDD-Find link in space gallery should be absent", () => shouldBe(page, el("BDD-Find link in space gallery"), "absent"));
      await session.step(46, "And BDD-Miss link in space gallery should be absent", () => shouldBe(page, el("BDD-Miss link in space gallery"), "absent"));
    });
    await run.scenario("Clearing the search brings both back", async () => {
      await session.step(49, "When user clears space search", () => clearField(page, el("space search")));
      await session.step(50, "Then BDD-Find link in space gallery should be visible", () => shouldBe(page, el("BDD-Find link in space gallery"), "visible"));
      await session.step(51, "And BDD-Miss link in space gallery should be visible", () => shouldBe(page, el("BDD-Miss link in space gallery"), "visible"));
    });
    await run.scenario("A child space is searchable inside its parent", async () => {
      await session.step(54, "When user picks \"Create Child Space...\" from the context menu of BDD-Find tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Find tree node inside browse tree")));
      await session.step(55, "And user enters \"BDD-Find-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Find-Child", el("Name input in Create Space dialog")));
      await session.step(56, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(57, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(58, "When user double-clicks on BDD-Find tree node inside browse tree", () => doubleClickOn(page, el("BDD-Find tree node inside browse tree")));
      await session.step(59, "Then BDD-Find-Child link in space gallery should be visible", () => shouldBe(page, el("BDD-Find-Child link in space gallery"), "visible"));
      await session.step(60, "When user enters \"zzz-no-such-space\" into space search", () => enterInto(page, "zzz-no-such-space", el("space search")));
      await session.step(61, "Then BDD-Find-Child link in space gallery should be absent", () => shouldBe(page, el("BDD-Find-Child link in space gallery"), "absent"));
      await session.step(62, "When user clears space search", () => clearField(page, el("space search")));
      await session.step(63, "Then BDD-Find-Child link in space gallery should be visible", () => shouldBe(page, el("BDD-Find-Child link in space gallery"), "visible"));
    });
    run.finish();
  });
});
