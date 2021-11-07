import * as ui from 'datagrok-api/ui'

import { renderBlock } from '../utils'
export default editor => {
    const blockManager = editor.BlockManager

    blockManager.add("basic-text-block", {
        category: "Datagrok",
        label: "Basic Text",
        icon: `<svg viewBox="0 0 50 40" fill="none" xmlns="http://www.w3.org/2000/svg">
          <path d="M23.0344 29L20.3544 21.88H8.35438L5.71438 29H0.634375L11.6744 0.639998H17.2344L28.2744 29H23.0344ZM14.3544 5.84L10.0344 17.44H18.6744L14.3544 5.84ZM30.2766 23.72C30.2766 22.0667 30.8099 20.7467 31.8766 19.76C32.9699 18.7467 34.3966 18.1067 36.1566 17.84L41.1566 17.08C42.1166 16.9467 42.5966 16.48 42.5966 15.68C42.5966 14.8 42.2899 14.0933 41.6766 13.56C41.0632 13 40.1432 12.72 38.9166 12.72C37.7432 12.72 36.8099 13.0533 36.1166 13.72C35.4232 14.36 35.0232 15.2 34.9166 16.24L30.6766 15.28C30.8632 13.52 31.7032 12.0267 33.1966 10.8C34.6899 9.54667 36.5832 8.92 38.8766 8.92C41.7032 8.92 43.7832 9.6 45.1166 10.96C46.4766 12.32 47.1566 14.0667 47.1566 16.2V25.88C47.1566 27.16 47.2366 28.2 47.3966 29H43.0766C42.9432 28.6 42.8766 27.7733 42.8766 26.52C41.5966 28.5733 39.5966 29.6 36.8766 29.6C34.9032 29.6 33.3032 29.0267 32.0766 27.88C30.8766 26.7067 30.2766 25.32 30.2766 23.72ZM37.7566 25.96C39.1966 25.96 40.3566 25.5733 41.2366 24.8C42.1432 24 42.5966 22.7067 42.5966 20.92V20.04L37.5166 20.8C35.7832 21.0933 34.9166 21.9733 34.9166 23.44C34.9166 24.1333 35.1699 24.7333 35.6766 25.24C36.1832 25.72 36.8766 25.96 37.7566 25.96Z" fill="#9A9797" />
        </svg>`,
        attributes: { class: "" },
        content: () => ui.divV([
            ui.h2('Insert title here'),
            ui.h3('Insert subtitle here'),
            ui.p('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua')
        ]).outerHTML,
        render: renderBlock,
        useBaseStyle: true
    })
}
