/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sticky-meta/schema-and-type.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.entity-type, views.entity-property-schemas]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, enterInto, expand, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountLower, rememberGalleryCount, stickyMetaFixturesGone, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("An entity type and a metadata schema, from creation to deletion", () => {
  const session = feature(test, "features/sticky-meta/schema-and-type.feature", import.meta.url);
  test("An entity type and a metadata schema, from creation to deletion", {tag: ["@journey", "@serial", "@sticky-meta", "@realizes:views.entity-type", "@realizes:views.entity-property-schemas"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the Sticky Meta schema \"bdd-sm-schema-{time}\" and entity type \"bdd-sm-type-{time}\" are removed now and at feature end", () => stickyMetaFixturesGone(page, session.text("bdd-sm-schema-{time}"), session.text("bdd-sm-type-{time}")));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(25, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(26, "And user expands \"Platform > Sticky Meta\" tree node inside browse tree", () => expand(page, el("\"Platform > Sticky Meta\" tree node inside browse tree")));
    await run.scenario("A new entity type needs a name and a matching expression (1.1)", async () => {
      await session.step(29, "When user clicks on \"Platform > Sticky Meta > Types\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Sticky Meta > Types\" tree node inside browse tree")));
      await session.step(30, "Then the \"Entity types\" view should be current", () => viewIsCurrent(page, "Entity types"));
      await session.step(31, "When user clicks on \"New Entity Type...\" button", () => clickOn(page, el("\"New Entity Type...\" button")));
      await session.step(32, "Then \"Create a new entity type\" dialog should be visible", () => shouldBe(page, el("\"Create a new entity type\" dialog"), "visible"));
      await session.step(33, "And OK button in \"Create a new entity type\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create a new entity type\" dialog"), "disabled"));
      await session.step(34, "When user enters \"bdd-sm-type-{time}\" into \"Name\" input in \"Create a new entity type\" dialog", () => enterInto(page, session.text("bdd-sm-type-{time}"), el("\"Name\" input in \"Create a new entity type\" dialog")));
      await session.step(35, "Then OK button in \"Create a new entity type\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create a new entity type\" dialog"), "disabled"));
      await session.step(36, "When user enters \"semtype=molecule\" into \"Matching expression\" input in \"Create a new entity type\" dialog", () => enterInto(page, "semtype=molecule", el("\"Matching expression\" input in \"Create a new entity type\" dialog")));
      await session.step(37, "Then OK button in \"Create a new entity type\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create a new entity type\" dialog"), "enabled"));
      await session.step(38, "When user clicks on OK button in \"Create a new entity type\" dialog", () => clickOn(page, el("OK button in \"Create a new entity type\" dialog")));
      await session.step(39, "Then the \"Create a new entity type\" dialog should close", () => dialogCloses(page, "Create a new entity type"));
      await session.step(40, "When user types \"bdd-sm-type-{time}\" into gallery search", () => typeInto(page, session.text("bdd-sm-type-{time}"), el("gallery search")));
      await session.step(41, "Then \"bdd-sm-type-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bdd-sm-type-{time}\" link in gallery")), "visible"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A schema is associated with the type and given four typed properties (1.2)", async () => {
      await session.step(46, "When user clicks on \"Platform > Sticky Meta > Schemas\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Sticky Meta > Schemas\" tree node inside browse tree")));
      await session.step(47, "Then the \"Schemas\" view should be current", () => viewIsCurrent(page, "Schemas"));
      await session.step(48, "When user clicks on \"New Schema...\" button", () => clickOn(page, el("\"New Schema...\" button")));
      await session.step(49, "Then \"Create a new schema\" dialog should be visible", () => shouldBe(page, el("\"Create a new schema\" dialog"), "visible"));
      await session.step(50, "And OK button in \"Create a new schema\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create a new schema\" dialog"), "disabled"));
      await session.step(51, "And \"Property Type\" input in \"Create a new schema\" dialog should offer \"string, int, bool, double, datetime, string_list\"", () => shouldOffer(page, el("\"Property Type\" input in \"Create a new schema\" dialog"), "string, int, bool, double, datetime, string_list"));
      await session.step(52, "When user enters \"bdd-sm-schema-{time}\" into \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, session.text("bdd-sm-schema-{time}"), el("\"Name\" input in \"Create a new schema\" dialog")));
      await session.step(53, "And user clicks on \"select entities\" action in \"Create a new schema\" dialog", () => clickOn(page, el("\"select entities\" action in \"Create a new schema\" dialog")));
      await session.step(54, "Then \"Select types for bdd-sm-schema-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Select types for bdd-sm-schema-{time}\" dialog")), "visible"));
      await session.step(55, "When user checks \"bdd-sm-type-{time}\" property in \"Select types for bdd-sm-schema-{time}\" dialog", () => check(page, el(session.text("\"bdd-sm-type-{time}\" property in \"Select types for bdd-sm-schema-{time}\" dialog"))));
      await session.step(56, "And user clicks on OK button in \"Select types for bdd-sm-schema-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Select types for bdd-sm-schema-{time}\" dialog"))));
      await session.step(57, "Then the \"Select types for bdd-sm-schema-{time}\" dialog should close", () => dialogCloses(page, session.text("Select types for bdd-sm-schema-{time}")));
      await session.step(58, "And \"Associated with\" input in \"Create a new schema\" dialog should contain text \"bdd-sm-type-{time}\"", () => shouldContainText(page, el("\"Associated with\" input in \"Create a new schema\" dialog"), session.text("bdd-sm-type-{time}")));
      await session.step(59, "When user enters \"rating\" into second \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "rating", el("second \"Name\" input in \"Create a new schema\" dialog")));
      await session.step(60, "And user selects \"int\" in \"Property Type\" input in \"Create a new schema\" dialog", () => selectIn(page, "int", el("\"Property Type\" input in \"Create a new schema\" dialog")));
      await session.step(61, "Then \"Property Type\" input in \"Create a new schema\" dialog should have value \"int\"", () => shouldHaveValue(page, el("\"Property Type\" input in \"Create a new schema\" dialog"), "int"));
      await session.step(62, "When user clicks on \"Add new property to schema\" button in \"Create a new schema\" dialog", () => clickOn(page, el("\"Add new property to schema\" button in \"Create a new schema\" dialog")));
      await session.step(63, "And user enters \"notes\" into third \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "notes", el("third \"Name\" input in \"Create a new schema\" dialog")));
      await session.step(64, "And user selects \"string\" in second \"Property Type\" input in \"Create a new schema\" dialog", () => selectIn(page, "string", el("second \"Property Type\" input in \"Create a new schema\" dialog")));
      await session.step(65, "Then second \"Property Type\" input in \"Create a new schema\" dialog should have value \"string\"", () => shouldHaveValue(page, el("second \"Property Type\" input in \"Create a new schema\" dialog"), "string"));
      await session.step(66, "When user clicks on \"Add new property to schema\" button in \"Create a new schema\" dialog", () => clickOn(page, el("\"Add new property to schema\" button in \"Create a new schema\" dialog")));
      await session.step(67, "And user enters \"verified\" into fourth \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "verified", el("fourth \"Name\" input in \"Create a new schema\" dialog")));
      await session.step(68, "And user selects \"bool\" in third \"Property Type\" input in \"Create a new schema\" dialog", () => selectIn(page, "bool", el("third \"Property Type\" input in \"Create a new schema\" dialog")));
      await session.step(69, "Then third \"Property Type\" input in \"Create a new schema\" dialog should have value \"bool\"", () => shouldHaveValue(page, el("third \"Property Type\" input in \"Create a new schema\" dialog"), "bool"));
      await session.step(70, "When user clicks on \"Add new property to schema\" button in \"Create a new schema\" dialog", () => clickOn(page, el("\"Add new property to schema\" button in \"Create a new schema\" dialog")));
      await session.step(71, "And user enters \"review_date\" into fifth \"Name\" input in \"Create a new schema\" dialog", () => enterInto(page, "review_date", el("fifth \"Name\" input in \"Create a new schema\" dialog")));
      await session.step(72, "And user selects \"datetime\" in fourth \"Property Type\" input in \"Create a new schema\" dialog", () => selectIn(page, "datetime", el("fourth \"Property Type\" input in \"Create a new schema\" dialog")));
      await session.step(73, "Then fourth \"Property Type\" input in \"Create a new schema\" dialog should have value \"datetime\"", () => shouldHaveValue(page, el("fourth \"Property Type\" input in \"Create a new schema\" dialog"), "datetime"));
      await session.step(74, "Then OK button in \"Create a new schema\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create a new schema\" dialog"), "enabled"));
      await session.step(75, "When user clicks on OK button in \"Create a new schema\" dialog", () => clickOn(page, el("OK button in \"Create a new schema\" dialog")));
      await session.step(76, "Then the \"Create a new schema\" dialog should close", () => dialogCloses(page, "Create a new schema"));
      await session.step(77, "When user types \"bdd-sm-schema-{time}\" into gallery search", () => typeInto(page, session.text("bdd-sm-schema-{time}"), el("gallery search")));
      await session.step(78, "Then \"bdd-sm-schema-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bdd-sm-schema-{time}\" link in gallery")), "visible"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
      await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Edit dialog shows the schema as it was saved (1.3)", async () => {
      await session.step(83, "When user picks \"Edit\" from the context menu of \"bdd-sm-schema-{time}\" link in gallery", () => pickFromContextMenu(page, "Edit", el(session.text("\"bdd-sm-schema-{time}\" link in gallery"))));
      await session.step(84, "Then \"Edit schema\" dialog should be visible", () => shouldBe(page, el("\"Edit schema\" dialog"), "visible"));
      await session.step(85, "And first \"Name\" input in \"Edit schema\" dialog should have value \"bdd-sm-schema-{time}\"", () => shouldHaveValue(page, el("first \"Name\" input in \"Edit schema\" dialog"), session.text("bdd-sm-schema-{time}")));
      await session.step(86, "And \"Associated with\" input in \"Edit schema\" dialog should contain text \"bdd-sm-type-{time}\"", () => shouldContainText(page, el("\"Associated with\" input in \"Edit schema\" dialog"), session.text("bdd-sm-type-{time}")));
      await session.step(87, "And second \"Name\" input in \"Edit schema\" dialog should have value \"rating\"", () => shouldHaveValue(page, el("second \"Name\" input in \"Edit schema\" dialog"), "rating"));
      await session.step(88, "And third \"Name\" input in \"Edit schema\" dialog should have value \"notes\"", () => shouldHaveValue(page, el("third \"Name\" input in \"Edit schema\" dialog"), "notes"));
      await session.step(89, "And fourth \"Name\" input in \"Edit schema\" dialog should have value \"verified\"", () => shouldHaveValue(page, el("fourth \"Name\" input in \"Edit schema\" dialog"), "verified"));
      await session.step(90, "And fifth \"Name\" input in \"Edit schema\" dialog should have value \"review_date\"", () => shouldHaveValue(page, el("fifth \"Name\" input in \"Edit schema\" dialog"), "review_date"));
      await session.step(91, "And sixth \"Name\" input in \"Edit schema\" dialog should be absent", () => shouldBe(page, el("sixth \"Name\" input in \"Edit schema\" dialog"), "absent"));
      await session.step(92, "When user clicks on CANCEL button in \"Edit schema\" dialog", () => clickOn(page, el("CANCEL button in \"Edit schema\" dialog")));
      await session.step(93, "Then the \"Edit schema\" dialog should close", () => dialogCloses(page, "Edit schema"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
      await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The schema is deleted from its list (1.4)", async () => {
      await session.step(98, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(99, "And user picks \"Delete\" from the context menu of \"bdd-sm-schema-{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"bdd-sm-schema-{time}\" link in gallery"))));
      await session.step(100, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(101, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(102, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(103, "And the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(104, "And \"bdd-sm-schema-{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"bdd-sm-schema-{time}\" link in gallery")), "absent"));
      await session.step(105, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(106, "Then no errors should have been logged", () => noErrors(page));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The entity type is deleted from its list (1.4)", async () => {
      await session.step(110, "When user clicks on \"Platform > Sticky Meta > Types\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Sticky Meta > Types\" tree node inside browse tree")));
      await session.step(111, "Then the \"Entity types\" view should be current", () => viewIsCurrent(page, "Entity types"));
      await session.step(112, "When user types \"bdd-sm-type-{time}\" into gallery search", () => typeInto(page, session.text("bdd-sm-type-{time}"), el("gallery search")));
      await session.step(113, "Then \"bdd-sm-type-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bdd-sm-type-{time}\" link in gallery")), "visible"));
      await session.step(114, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(115, "And user picks \"Delete\" from the context menu of \"bdd-sm-type-{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"bdd-sm-type-{time}\" link in gallery"))));
      await session.step(116, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(117, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(118, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(119, "And the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(120, "And \"bdd-sm-type-{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"bdd-sm-type-{time}\" link in gallery")), "absent"));
      await session.step(121, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(122, "Then no errors should have been logged", () => noErrors(page));
      await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
