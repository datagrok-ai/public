import * as ui from 'datagrok-api/ui'

import { renderBlock } from '../utils'
export default editor => {
  const blockManager = editor.BlockManager

  blockManager.add("headline-block", {
    category: "Datagrok",
    label: "Headline",
    icon: `<svg viewBox="10 0 47.5 40" fill="none"  xmlns="http://www.w3.org/2000/svg">
            <path d="M46.9772 29H42.1372V16.88H28.8172V29H24.0172V0.639998H28.8172V12.36H42.1372V0.639998H46.9772V29Z" fill="#9A9797" />
          </svg>`,
    attributes: {},
    content: () => ui.h1('This is a headline text').outerHTML,
    render: renderBlock,
    useBaseStyle: true
  })
}
