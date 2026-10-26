/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/gallery-search-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Gallery search edge cases: characters of the search syntax typed as text", () => {
  const session = feature(test, "features/browse/gallery-search-edge-cases.feature", import.meta.url);
  test("A search of an opening parenthesis empties the list with no error, and clearing it brings the space back [what=an opening parenthesis, text=(]", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And a space named \"BDDEdgeSearch{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSearch{time}")));
    await session.step(17, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(18, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
    await session.step(19, "And user enters \"BDDEdgeSearch{time}\" into gallery search", () => enterInto(page, session.text("BDDEdgeSearch{time}"), el("gallery search")));
    await session.step(20, "Then BDDEdgeSearch{time} link in gallery should be visible", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "visible"));
    await session.step(23, "When user enters \"(\" into gallery search", () => enterInto(page, "(", el("gallery search")));
    await session.step(24, "Then BDDEdgeSearch{time} link in gallery should be absent", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "absent"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(27, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(28, "Then BDDEdgeSearch{time} link in gallery should be visible", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "visible"));
  });
  test("A search of a double quote empties the list with no error, and clearing it brings the space back [what=a double quote, text=\\\"]", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And a space named \"BDDEdgeSearch{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSearch{time}")));
    await session.step(17, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(18, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
    await session.step(19, "And user enters \"BDDEdgeSearch{time}\" into gallery search", () => enterInto(page, session.text("BDDEdgeSearch{time}"), el("gallery search")));
    await session.step(20, "Then BDDEdgeSearch{time} link in gallery should be visible", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "visible"));
    await session.step(23, "When user enters \"\\\"\" into gallery search", () => enterInto(page, "\"", el("gallery search")));
    await session.step(24, "Then BDDEdgeSearch{time} link in gallery should be absent", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "absent"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(27, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(28, "Then BDDEdgeSearch{time} link in gallery should be visible", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "visible"));
  });
  test("A search of a condition with no value empties the list with no error, and clearing it brings the space back [what=a condition with no value, text=name =]", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And a space named \"BDDEdgeSearch{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSearch{time}")));
    await session.step(17, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(18, "When user clicks on Spaces tree node inside browse tree", () => clickOn(page, el("Spaces tree node inside browse tree")));
    await session.step(19, "And user enters \"BDDEdgeSearch{time}\" into gallery search", () => enterInto(page, session.text("BDDEdgeSearch{time}"), el("gallery search")));
    await session.step(20, "Then BDDEdgeSearch{time} link in gallery should be visible", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "visible"));
    await session.step(23, "When user enters \"name =\" into gallery search", () => enterInto(page, "name =", el("gallery search")));
    await session.step(24, "Then BDDEdgeSearch{time} link in gallery should be absent", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "absent"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(27, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(28, "Then BDDEdgeSearch{time} link in gallery should be visible", () => shouldBe(page, el(session.text("BDDEdgeSearch{time} link in gallery")), "visible"));
  });
});
