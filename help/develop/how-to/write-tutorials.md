<!-- TITLE: Write tutorials -->

# Tutorials

## Definition

A *tutorial* is a collection of interactive materials that educate users about
the platform. Such set of instructions can explain key concepts, such as
[viewers](../../visualize/viewers.md) and
[filters](../../visualize/viewers/filters.md), or illustrate certain workflows,
e.g., how to [access data](../../access/access.md) or train a [predictive
model](../../learn/predictive-modeling.md). A collection of tutorials on a
specific topic is called a *track*.

## Writing a tutorial

From a programmer's point of view, a tutorial is a class containing certain
attributes as well as an event-based set of instructions. To start writing a
tutorial in a package, we need to install `@datagrok-libraries/tutorials`.

The first step is to declare a subclass of `Tutorial`:

```typescript
class CustomTutorial extends Tutorial {
  // ...
}
```

Each tutorial should have `name` and `description` as well as an implementation
of the `_run` method, where you should define a number of steps. Steps are the
key units in every tutorial, and they correspond to the `action` method or
helper methods calling it internally. Steps can come with interactive hints or
more detailed description.

## Registering a tutorial

Tutorials from outside the public
[Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials)
package should be registered. To do that, add a function to your `package.ts`
file. It should have the `tutorial` tag and the `meta.name` parameter. As it
should return a tutorial instance, the output type must be `object`. There are
also optional parameters that you may use: `meta.track` specifies the track
the tutorial belongs to; `meta.icon` contains the icon path relative to the
package root; `description` contents are displayed in the UI.

```typescript
import { CustomTutorial } from './tutorials/custom-tutorial';


//tags: tutorial
//meta.icon: images/custom-tutorial.png
//meta.name: Custom Tutorial
//meta.track: Test Track
//description: This tutorial illustrates a new feature
//output: object
export function registerTutorial() {
  return new CustomTutorial();
}
```

Without a specified track, a tutorial is assigned to a new track created
specially for the package. If an existing track is mentioned in the annotation,
a tutorial is assigned to it. So, it's possible to extend a track from another
package. If the `meta.icon` tag is absent, the default icon is used. If there is
no description for the tutorial, it gets added with empty description.

It is also possible to register a track:

```typescript
//tags: track
//help-url: https://datagrok.ai/help
//output: object
//meta.name: Test Track
export function registerTrack() {
  return new Track('Test Track', [], 'https://datagrok.ai/help');
}
```

The function should have the `track` tag, parameters `help-url` and `meta.name`,
and it should return a `Track` instance (marked as `output: object` in the
function header).

## Standard tutorials visibility

The `Tutorials` package provides a set of default tracks, such as *Exploratory
Data Analysis* or *Data Access*. You can control their visibility on the user
group level by applying relevant [package settings](../develop.md#package-settings).

See also:

* [Packages](../develop.md#packages)
* [JavaScript API](../js-api.md)
