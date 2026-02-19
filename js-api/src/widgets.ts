/**
 * @module widgets
 *
 * UI widgets and components for building Datagrok applications.
 *
 * ## Module Organization
 *
 * This module contains 50+ exports organized into logical categories:
 *
 * ### Core Widget Infrastructure
 * - {@link Widget} - Base class for all JavaScript-based widgets
 * - {@link DartWidget} - Base class for widgets backed by Dart implementation
 * - {@link ObjectPropertyBag} - Proxy-based property accessor for dynamic property access
 * - {@link DartWrapper} - Low-level Dart object wrapper
 *
 * ### Container & Layout Widgets
 * - {@link Accordion}, {@link AccordionPane} - Collapsible panel container
 * - {@link TabControl}, {@link TabPane} - Tabbed interface container
 * - {@link ToolboxPage} - Toolbox page for organizing tools
 *
 * ### Menu System
 * - {@link Menu} - Context menus and menu bars with rich options
 * - {@link Balloon} - Tooltip/popup balloon notifications
 * - Menu interfaces: {@link IMenuItemsOptions}, {@link IMenuColorPaletteOptions},
 *   {@link IMenuFontEditorOptions}, {@link IMenuColumnSelectorOptions}, etc.
 *
 * ### Input Controls
 * - {@link InputBase} - Base class for all input controls
 * - {@link InputForm} - Form container for multiple inputs
 * - {@link DateInput} - Date/time picker input
 * - {@link ChoiceInput} - Dropdown/select input
 * - {@link TagEditor}, {@link TagElement} - Tag/chip input editor
 * - {@link TypeAhead} - Autocomplete text input
 * - {@link MarkdownInput} - Markdown text editor
 * - {@link CodeInput}, {@link CodeEditor} - Code editing with syntax highlighting
 *
 * ### Tree & List Components
 * - {@link TreeViewNode}, {@link TreeViewGroup} - Hierarchical tree view
 * - {@link HtmlTable} - HTML table widget
 *
 * ### Data Visualization & Selection
 * - {@link RangeSlider} - Range selection slider
 * - {@link ColumnComboBox} - Column selector dropdown
 * - {@link Legend} - Chart legend widget
 * - {@link PropertyGrid} - Property editor grid
 *
 * ### Specialized Widgets
 * - {@link FilesWidget} - File browser widget
 * - {@link Breadcrumbs} - Navigation breadcrumbs
 * - {@link DropDown} - Generic dropdown component
 * - {@link FunctionsWidget} - Functions browser widget
 * - {@link Favorites} - Favorites management widget
 * - {@link VisualDbQueryEditor} - Visual database query editor
 * - {@link Dialog} - Modal dialog window
 *
 * ### Progress Indicators
 * - {@link TaskBarProgressIndicator} - Progress indicator in task bar
 *
 * ## Usage Patterns
 *
 * Most widgets follow these patterns:
 *
 * 1. **Dart-backed widgets** extend {@link DartWidget} and wrap a Dart object:
 *    ```typescript
 *    class MyWidget extends DartWidget {
 *      constructor(dart: any) { super(dart); }
 *    }
 *    ```
 *
 * 2. **Property access** via {@link ObjectPropertyBag} proxy:
 *    ```typescript
 *    viewer.props.showAxes = true;  // Dynamic property access
 *    ```
 *
 * 3. **Event subscriptions** via RxJS Observables:
 *    ```typescript
 *    widget.onChanged.subscribe(e => console.log('Changed:', e));
 *    ```
 *
 * @see {@link ../ui.ts} for UI factory functions that create these widgets
 */

export * from './widgets/index';
