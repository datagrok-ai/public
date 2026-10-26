/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sticky-meta/schema-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.entity-property-schemas, views.entity-type, core#677]
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
import {check, clickOn, enterInto, expand, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, contextPanelOpen, dialogCloses, noEntityType, noStickyMetaSchema, stickyMetaFixturesGone, stickyMetaSchemaHeld, stickyMetaSchemaOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, infoBalloon, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Sticky Meta schema edge cases: a taken name, a malformed matching expression, a bad address, delete from the tree", () => {
  const session = feature(test, "features/sticky-meta/schema-edge-cases.feature", import.meta.url);
  test("A new schema under a taken name asks before dropping its properties, and CANCEL keeps them", {tag: ["@serial", "@sticky-meta", "@edge-cases", "@realizes:views.entity-property-schemas", "@realizes:views.entity-type"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And a Sticky Meta schema \"bdd-sme-schema-{time}\" for molecules, with the entity type \"bdd-sme-type-{time}\", is on the server with the properties:", () => stickyMetaSchemaOnServer(page, session.text("bdd-sme-schema-{time}"), session.text("bdd-sme-type-{time}"), [["rating","int"],["notes","string"]]), [["rating","int"],["notes","string"]]);
    await session.step(23, "And the Sticky Meta schema \"bdd-sme-none-{time}\" and entity type \"bdd-sme-badtype-{time}\" are removed now and at feature end", () => stickyMetaFixturesGone(page, session.text("bdd-sme-none-{time}"), session.text("bdd-sme-badtype-{time}")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(27, "And user expands \"Platform > Sticky Meta\" tree node inside browse tree", () => expand(page, el("\"Platform > Sticky Meta\" tree node inside browse tree")));
    await session.step(30, "When user clicks on \"Platform > Sticky Meta > Schemas\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Sticky Meta > Schemas\" tree node inside browse tree")));
    await session.step(31, "Then the \"Schemas\" view should be current", () => viewIsCurrent(page, "Schemas"));
    await session.step(32, "When user clicks on \"New Schema...\" button", () => clickOn(page, el("\"New Schema...\" button")));
    await session.step(33, "Then \"Create a new schema\" dialog should be visible", () => shouldBe(page, el("\"Create a new schema\" dialog"), "visible"));
    await session.step(34, "When user enters \"bdd-sme-schema-{time}\" into \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, session.text("bdd-sme-schema-{time}"), el("\"Name\" input in \"Create a new schema\" dialog")));
    await session.step(35, "And user clicks on \"select entities\" action in \"Create a new schema\" dialog", () => clickOn(page, el("\"select entities\" action in \"Create a new schema\" dialog")));
    await session.step(36, "Then \"Select types for bdd-sme-schema-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Select types for bdd-sme-schema-{time}\" dialog")), "visible"));
    await session.step(37, "When user checks \"bdd-sme-type-{time}\" property in \"Select types for bdd-sme-schema-{time}\" dialog", () => check(page, el(session.text("\"bdd-sme-type-{time}\" property in \"Select types for bdd-sme-schema-{time}\" dialog"))));
    await session.step(38, "And user clicks on OK button in \"Select types for bdd-sme-schema-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Select types for bdd-sme-schema-{time}\" dialog"))));
    await session.step(39, "Then the \"Select types for bdd-sme-schema-{time}\" dialog should close", () => dialogCloses(page, session.text("Select types for bdd-sme-schema-{time}")));
    await session.step(40, "When user enters \"score\" into second \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "score", el("second \"Name\" input in \"Create a new schema\" dialog")));
    await session.step(41, "And user selects \"int\" in \"Property Type\" input in \"Create a new schema\" dialog", () => selectIn(page, "int", el("\"Property Type\" input in \"Create a new schema\" dialog")));
    await session.step(42, "Then OK button in \"Create a new schema\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create a new schema\" dialog"), "enabled"));
    await session.step(43, "When user clicks on OK button in \"Create a new schema\" dialog", () => clickOn(page, el("OK button in \"Create a new schema\" dialog")));
    await session.step(44, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(45, "And \"Are you sure?\" dialog should contain text \"delete properties rating, notes\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "delete properties rating, notes"));
    await session.step(46, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
    await session.step(47, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(48, "And \"Create a new schema\" dialog should be visible", () => shouldBe(page, el("\"Create a new schema\" dialog"), "visible"));
    await session.step(49, "And the Sticky Meta schema \"bdd-sme-schema-{time}\" should be on the server with the properties \"rating: int, notes: string\"", () => stickyMetaSchemaHeld(page, session.text("bdd-sme-schema-{time}"), "rating: int, notes: string"));
    await session.step(50, "When user clicks on CANCEL button in \"Create a new schema\" dialog", () => clickOn(page, el("CANCEL button in \"Create a new schema\" dialog")));
    await session.step(51, "Then the \"Create a new schema\" dialog should close", () => dialogCloses(page, "Create a new schema"));
    await session.step(52, "And the Sticky Meta schema \"bdd-sme-schema-{time}\" should be on the server with the properties \"rating: int, notes: string\"", () => stickyMetaSchemaHeld(page, session.text("bdd-sme-schema-{time}"), "rating: int, notes: string"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A matching expression not in tag=value form is refused", {tag: ["@serial", "@sticky-meta", "@edge-cases", "@realizes:views.entity-property-schemas", "@realizes:views.entity-type"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And a Sticky Meta schema \"bdd-sme-schema-{time}\" for molecules, with the entity type \"bdd-sme-type-{time}\", is on the server with the properties:", () => stickyMetaSchemaOnServer(page, session.text("bdd-sme-schema-{time}"), session.text("bdd-sme-type-{time}"), [["rating","int"],["notes","string"]]), [["rating","int"],["notes","string"]]);
    await session.step(23, "And the Sticky Meta schema \"bdd-sme-none-{time}\" and entity type \"bdd-sme-badtype-{time}\" are removed now and at feature end", () => stickyMetaFixturesGone(page, session.text("bdd-sme-none-{time}"), session.text("bdd-sme-badtype-{time}")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(27, "And user expands \"Platform > Sticky Meta\" tree node inside browse tree", () => expand(page, el("\"Platform > Sticky Meta\" tree node inside browse tree")));
    await session.step(57, "When user clicks on \"Platform > Sticky Meta > Types\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Sticky Meta > Types\" tree node inside browse tree")));
    await session.step(58, "Then the \"Entity types\" view should be current", () => viewIsCurrent(page, "Entity types"));
    await session.step(59, "When user clicks on \"New Entity Type...\" button", () => clickOn(page, el("\"New Entity Type...\" button")));
    await session.step(60, "Then \"Create a new entity type\" dialog should be visible", () => shouldBe(page, el("\"Create a new entity type\" dialog"), "visible"));
    await session.step(61, "When user enters \"bdd-sme-badtype-{time}\" into \"Name\" input in \"Create a new entity type\" dialog", () => enterInto(page, session.text("bdd-sme-badtype-{time}"), el("\"Name\" input in \"Create a new entity type\" dialog")));
    await session.step(62, "And user enters \"semtype\" into \"Matching expression\" input in \"Create a new entity type\" dialog", () => enterInto(page, "semtype", el("\"Matching expression\" input in \"Create a new entity type\" dialog")));
    await session.step(63, "Then \"Matching expression\" input in \"Create a new entity type\" dialog should be invalid", () => shouldBe(page, el("\"Matching expression\" input in \"Create a new entity type\" dialog"), "invalid"));
    await session.step(64, "And OK button in \"Create a new entity type\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create a new entity type\" dialog"), "disabled"));
    await session.step(65, "When user enters \"semtype=\" into \"Matching expression\" input in \"Create a new entity type\" dialog", () => enterInto(page, "semtype=", el("\"Matching expression\" input in \"Create a new entity type\" dialog")));
    await session.step(66, "Then \"Matching expression\" input in \"Create a new entity type\" dialog should be invalid", () => shouldBe(page, el("\"Matching expression\" input in \"Create a new entity type\" dialog"), "invalid"));
    await session.step(67, "And OK button in \"Create a new entity type\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create a new entity type\" dialog"), "disabled"));
    await session.step(68, "When user enters \"semtype=molecule\" into \"Matching expression\" input in \"Create a new entity type\" dialog", () => enterInto(page, "semtype=molecule", el("\"Matching expression\" input in \"Create a new entity type\" dialog")));
    await session.step(69, "Then \"Matching expression\" input in \"Create a new entity type\" dialog should be valid", () => shouldBe(page, el("\"Matching expression\" input in \"Create a new entity type\" dialog"), "valid"));
    await session.step(70, "And OK button in \"Create a new entity type\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create a new entity type\" dialog"), "enabled"));
    await session.step(71, "When user clicks on CANCEL button in \"Create a new entity type\" dialog", () => clickOn(page, el("CANCEL button in \"Create a new entity type\" dialog")));
    await session.step(72, "Then the \"Create a new entity type\" dialog should close", () => dialogCloses(page, "Create a new entity type"));
    await session.step(73, "And no entity type named \"bdd-sme-badtype-{time}\" should be on the server", () => noEntityType(page, session.text("bdd-sme-badtype-{time}")));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
    await session.step(75, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An address naming no schema says so and leaves the current view", {tag: ["@serial", "@sticky-meta", "@edge-cases", "@realizes:views.entity-property-schemas", "@realizes:views.entity-type"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And a Sticky Meta schema \"bdd-sme-schema-{time}\" for molecules, with the entity type \"bdd-sme-type-{time}\", is on the server with the properties:", () => stickyMetaSchemaOnServer(page, session.text("bdd-sme-schema-{time}"), session.text("bdd-sme-type-{time}"), [["rating","int"],["notes","string"]]), [["rating","int"],["notes","string"]]);
    await session.step(23, "And the Sticky Meta schema \"bdd-sme-none-{time}\" and entity type \"bdd-sme-badtype-{time}\" are removed now and at feature end", () => stickyMetaFixturesGone(page, session.text("bdd-sme-none-{time}"), session.text("bdd-sme-badtype-{time}")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(27, "And user expands \"Platform > Sticky Meta\" tree node inside browse tree", () => expand(page, el("\"Platform > Sticky Meta\" tree node inside browse tree")));
    await session.step(78, "When user clicks on \"Platform > Sticky Meta > Schemas\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Sticky Meta > Schemas\" tree node inside browse tree")));
    await session.step(79, "Then the \"Schemas\" view should be current", () => viewIsCurrent(page, "Schemas"));
    await session.step(80, "When user opens the address \"/meta/schemas/bdd-sme-none-{time}\"", () => openAddress(page, session.text("/meta/schemas/bdd-sme-none-{time}")));
    await session.step(81, "Then an error balloon containing \"Schema \\\"bdd-sme-none-{time}\\\" not found\" should have been shown", () => errorBalloonText(page, session.text("Schema \"bdd-sme-none-{time}\" not found")));
    await session.step(82, "And the \"Schemas\" view should be current", () => viewIsCurrent(page, "Schemas"));
  });
  test("A schema deleted from the Browse tree leaves the tree at once", {tag: ["@serial", "@sticky-meta", "@edge-cases", "@realizes:views.entity-property-schemas", "@realizes:views.entity-type", "@known-failure:95", "@realizes:core#677"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And a Sticky Meta schema \"bdd-sme-schema-{time}\" for molecules, with the entity type \"bdd-sme-type-{time}\", is on the server with the properties:", () => stickyMetaSchemaOnServer(page, session.text("bdd-sme-schema-{time}"), session.text("bdd-sme-type-{time}"), [["rating","int"],["notes","string"]]), [["rating","int"],["notes","string"]]);
    await session.step(23, "And the Sticky Meta schema \"bdd-sme-none-{time}\" and entity type \"bdd-sme-badtype-{time}\" are removed now and at feature end", () => stickyMetaFixturesGone(page, session.text("bdd-sme-none-{time}"), session.text("bdd-sme-badtype-{time}")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(27, "And user expands \"Platform > Sticky Meta\" tree node inside browse tree", () => expand(page, el("\"Platform > Sticky Meta\" tree node inside browse tree")));
    await knownFailure(async () => {
      await session.step(88, "When user expands \"Platform > Sticky Meta > Schemas\" tree node inside browse tree", () => expand(page, el("\"Platform > Sticky Meta > Schemas\" tree node inside browse tree")));
      await session.step(89, "Then \"Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}\" tree node inside browse tree")), "visible"));
      await session.step(90, "When user picks \"Delete\" from the context menu of \"Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}\" tree node inside browse tree", () => pickFromContextMenu(page, "Delete", el(session.text("\"Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}\" tree node inside browse tree"))));
      await session.step(91, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(92, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(93, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(94, "And no Sticky Meta schema named \"bdd-sme-schema-{time}\" should be on the server", () => noStickyMetaSchema(page, session.text("bdd-sme-schema-{time}")));
      await session.step(95, "And \"Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"Platform > Sticky Meta > Schemas > bdd-sme-schema-{time}\" tree node inside browse tree")), "absent"));
      await session.step(96, "And an info balloon should have been shown", () => infoBalloon(page));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    }, 95);
  });
});
