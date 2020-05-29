import {Viewer} from "./src/viewer";
import {VirtualView} from "./src/view";
import {Accordion, Dialog, InputBase, Menu, TabControl} from "./src/ui_classes";

export function e (s, cl = null) {
        let x = document.createElement(s);
        if (cl !== null)
            _class(x, cl);
        return x;
    }

export function appendAll(root, elements) {
        let fragment = document.createDocumentFragment();
        for (let i = 0; i < elements.length; i++) {
            let e = elements[i];
            if (e instanceof Viewer)
                e = e.root;
            fragment.appendChild(e);
        }
        root.appendChild(fragment);
        return root;
    }

export function _innerText(x, s) { x.innerText = s; return x; }
export function _class(x, s) { x.classList.add(s); return x; }
export function _color(x, s) { x.style.color = s; return x; }
export function _backColor(x, s) { x.style.backgroundColor = s; return x; }

export function canvas() { return document.createElement("CANVAS"); }

export function h1(s) { return _innerText(ui.e('h1'), s); }
export function h2(s) { let x = e('h2'); x.innerText = s; return x; }
export function h3(s) { let x = e('h3'); x.innerText = s; return x; }

export function accordion() { return Accordion.create(); }
export function tabControl(vertical = false) { return TabControl.create(vertical); }

export function divText(s) {
        let e = document.createElement('div');
        e.innerText = s;
        return e;
    }

export function iconFA(name, handler, tooltipMsg = null) {
        let i = document.createElement('i');
        i.classList.add('grok-icon');
        i.classList.add('fal');
        i.classList.add(`fa-${name}`);
        i.addEventListener("click", handler);
        if (tooltipMsg !== null)
            tooltip(i, tooltipMsg);
        return i;
    }

export function render(x) { return grok_UI_Render(x); }
export function renderCard(x) { return grok_UI_RenderCard(x); }
export function span(x) { return grok_UI_Span(x); }

export function div(items = [], className = null) { return grok_UI_Div(items, className); }
export function divV(items, className = null) { return grok_UI_DivV(items, className); }
export function divH(items, className = null) { return grok_UI_DivH(items, className); }

export function card(content) { return ui.div([content], 'd4-item-card'); }

export function loader() { return grok_UI_Loader(); }
export function setUpdateIndicator(element, updating = true) { return grok_UI_SetUpdateIndicator(element, updating)};

export function button(text, handler, tooltip = null) { return grok_UI_Button(text, handler, tooltip); }
export function bigButton(text, handler, tooltip = null) { return grok_UI_BigButton(text, handler, tooltip); }

/** Creates a visual table based on [map]. */
export function tableFromMap(map) { return grok_UI_TableFromMap(map); }

/** Creates a visual element representing list of [items]. */
export function list(items) { return grok_UI_List(Array.from(items).map(_toDart)); }

/** Creates a [Dialog].
 * @returns {Dialog} */
export function dialog(title = '') { return Dialog.create(title); }

/** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
 * tooltip, and popup menu.
 * Returns [element]. */
export function bind(item, element) { return grok_UI_Bind(item, element); }

export function virtualView(length, renderer) {
        let view = VirtualView.create();
        view.setData(length, renderer);
        return view.root;
    }

export function popupMenu(items) {
        function populate(menu, item) {
            for (let key of Object.keys(item)) {
                let value = item[key];
                if (value instanceof Function)
                    menu.item(key, value);
                else
                    populate(menu.group(key), value);
            }
        }

        let menu = Menu.popup();
        populate(menu, items);
        menu.show();
    }

export function tooltipHide() { grok_Tooltip_Hide(); }

export function tooltip(e, x) { grok_Tooltip_SetOn(e, x); return e; }

export function tooltipShow(content, x, y) { grok_Tooltip_Show(content, x, y); }

export function inputs(inputs) { return div(inputs.map((x) => x.root), 'pure-form,pure-form-aligned');}

export function intInput(name, value) { return new InputBase(grok_IntInput(name, value)); }
export function choiceInput(name, selected, items) { return new InputBase(grok_ChoiceInput(name, selected, items)); }
export function multiChoiceInput(name, value, items) { return new InputBase(grok_MultiChoiceInput(name, value, items)); }
export function stringInput(name, value) { return new InputBase(grok_StringInput(name, value)); }
export function floatInput(name, value) { return new InputBase(grok_FloatInput(name, value)); }
export function dateInput(name, value) { return new InputBase(grok_DateInput(name, value.d)); }
export function boolInput(name, value, callback = null) {
    return new InputBase(grok_BoolInput(name, value), callback);
}
export function moleculeInput(name, value) { return new InputBase(grok_MoleculeInput(name, value)); }
export function columnInput(name, table, value) { return new InputBase(grok_ColumnInput(name, table.d, value.d)); }
export function columnsInput(name, table) { return new InputBase(grok_ColumnsInput(name, table.d)); }



/** UI Tools **/
export class tools {

    static handleResize(element, onChanged) {
        let width = element.clientWidth;
        let height = element.clientHeight;
        let interval = setInterval(() => {
            let newWidth = element.clientWidth;
            let newHeight = element.clientHeight;
            if (newWidth !== width || newHeight !== height) {
                width = newWidth;
                height = newHeight;
                onChanged(width, height);
            }
        }, 100);
        return () => clearInterval(interval);
    }
}

