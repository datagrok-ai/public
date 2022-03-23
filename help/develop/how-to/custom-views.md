<!-- TITLE: Create a custom view -->

# Custom views

A view is a set of visualizations grouped together. Typically it is associated with a particular dataframe (then it is
called a [table view](../../overview/table-view.md)). However, essentially a view can contain pretty much anything,
that's what we will discuss in a while. When a view is registered on the platform, you can do the following:

* Open and pass parameters to the view from the URL
* Save the view as part of a project
* Add the view to the navigation bar
* Link the view with [entities](../../overview/objects.md) (custom or the default ones)

Let's start with the simplest example:

```javascript
let view = grok.shell.newView('Custom View', [ui.divText('Hi!')]);
```

In just one line of code we have created a view instance that you can later fill with new visual components, such as
icons, buttons, or side panels. See how it works (and looks!) in
the [full example](https://public.datagrok.ai/js/samples/ui/views/views).

As simple and easy as this method may seem, it is still more suitable for ad-hoc experiments. To implement a view from a
package, you would need to extend the [ViewBase](https://datagrok.ai/js-api/classes/dg.ViewBase) class. For
instance, the following part of code defines a new view for Jupyter Notebooks:

```javascript
class NotebookView extends DG.ViewBase {
    constructor(params, path) {
        super(params, path);
        this.TYPE = 'Notebook';
        this.PATH = '/notebook';
    }

    // Override basic methods
    get type() { return this.TYPE };
    get helpUrl() { return '/help/compute/jupyter-notebook.md'; }
    get name() { return 'Notebook' };
    get path() { return `${this.PATH}/${this.notebookId}` };

    // Icon
    getIcon() {
        let img = document.createElement('img');
        img.src = '/images/entities/jupyter.png';
        img.height = 18;
        img.width = 18;
        return img;
    };

    // View state serialization/deserialization
    saveStateMap() { return {'notebookId': this.notebookId }; }
    loadStateMap(stateMap) { open(stateMap['notebookId']); }

    // URL path handler
    handlePath(path) {
        let id = path.replace(`${this.PATH}/`, '');
        open(id);
    }

    // URL path checker
    acceptsPath(path) { return path.startsWith(this.PATH); }
}
```

The function `notebookView` defined in the main package file allocates an instance of `NotebookView`:

```javascript
//name: Notebook
//description: Creates a Notebook View
//input: map params
//input: string path
//tags: view
//output: view result
export function notebookView(params = null, path = '') {
    return new NotebookView(params, path);
}
```

Follow these function annotation rules to register a view:

* the header parameters should contain the `view` tag
* add two inputs `params` and `path` to pass URL parameters and path respectively
* specify an output of the `view` type.

These are the fragments of real code from the
public [Notebooks](https://github.com/datagrok-ai/public/blob/master/packages/Notebooks/src/package.js)
package. This is what allows users to work with Jupyter Notebooks on Datagrok (see examples in
the [Notebooks](https://public.datagrok.ai/notebooks?) browser).

See also:

* [JavaScript API](../js-api.md)
* [JavaScript API: View](https://datagrok.ai/js-api/classes/dg.View)
* [JavaScript API: ViewBase](https://datagrok.ai/js-api/classes/dg.ViewBase)
* [JavaScript API Samples: Custom view](https://public.datagrok.ai/js/samples/ui/views/views)
* [JavaScript API Samples: Virtual view](https://public.datagrok.ai/js/samples/ui/virtual-view)
* [Table view](../../overview/table-view.md)
* [View layout](../../visualize/view-layout.md)
* [Routing](../../overview/routing.md)
