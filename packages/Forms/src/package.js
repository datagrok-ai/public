/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok'
import * as ui from 'datagrok-api/ui'
import * as DG from 'datagrok-api/dg'

/* 3rd party imports */
import grapesjs from 'grapesjs'
import 'grapesjs-blocks-basic'
import 'grapesjs/dist/css/grapes.min.css'

/* Project imports */
import '../css/style.css'

/* constants */
export const _package = new DG.Package()

//name: FormsApp
//tags: app
export function FormsApp() {
  const view = grok.shell.newView('Form Builder')
  const rootEl = view.root
  rootEl.setAttribute('id', 'gjs-root')
  rootEl.style.backgroundColor = '#2a2a2a'

  const stylesContainer = ui.panel([], { id: 'styles-container', style: { backgroundColor: '#444', padding: '0' } })
  const layersContainer = ui.panel([], { id: 'layers-container', style: { backgroundColor: '#444', padding: '0' } })
  const blocksContainer = ui.panel([], { id: 'blocks-container', style: { backgroundColor: '#444', padding: '0' } })
  const traitsContainer = ui.panel([], { id: 'traits-container', style: { backgroundColor: '#444', padding: '0' } })

  const tabControl = ui.tabControl({
    'STYLES': stylesContainer,
    'LAYERS': layersContainer,
    'BLOCKS': blocksContainer,
    'TRAITS': traitsContainer
  })
  tabControl.root.style.height = '100%'

  const propertyPanel = document.getElementsByClassName('grok-prop-panel').item(0)
  propertyPanel.style.display = 'block'
  propertyPanel.innerHTML = ''
  propertyPanel.append(tabControl.root)

  const gjsEditor = grapesjs.init({
    container: '#gjs-root',
    plugins: ["gjs-blocks-basic"],
    height: '100%',
    width: '100%',
    storageManager: false,

    panels: {
      defaults: []
    },
    selectorManager: {
      appendTo: stylesContainer
    },

    layerManager: {
      appendTo: layersContainer
    },

    styleManager: {
      appendTo: stylesContainer,
      sectors: [{
        name: 'Dimension',
        open: false,
        buildProps: ['width', 'min-height', 'padding'],
        properties: [
          {
            type: 'integer',
            name: 'The width',
            property: 'width',
            units: ['px', '%'],
            defaults: 'auto',
            min: 0,
          }
        ]
      }, {
        name: 'Extra',
        open: false,
        buildProps: ['background-color', 'box-shadow', 'custom-prop'],
        properties: [
          {
            id: 'custom-prop',
            name: 'Custom Label',
            property: 'font-size',
            type: 'select',
            defaults: '32px',
            options: [
              { value: '12px', name: 'Tiny' },
              { value: '18px', name: 'Medium' },
              { value: '32px', name: 'Big' },
            ],
          }
        ]
      }]
    },
    blockManager: {
      appendTo: blocksContainer,
      blocks: [
        {
          id: 'section', // id is mandatory
          label: '<b>Section</b>', // You can use HTML/SVG inside labels
          attributes: { class: 'gjs-block-section' },
          content: `<section>
            <h1>This is a simple title</h1>
            <div>This is just a Lorem text: Lorem ipsum dolor sit amet</div>
          </section>`,
        }, {
          id: 'text',
          label: 'Text',
          content: '<div data-gjs-type="text">Insert your text here</div>',
        }, {
          id: 'image',
          label: 'Image',
          // Select the component once it's dropped
          select: true,
          // You can pass components as a JSON instead of a simple HTML string,
          // in this case we also use a defined component type `image`
          content: { type: 'image' },
          // This triggers `active` event on dropped components and the `image`
          // reacts by opening the AssetManager
          activate: true,
        }
      ]
    },
    traitManager: {
      appendTo: traitsContainer,
    },
  }
  )
}
