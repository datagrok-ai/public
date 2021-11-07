/**
 * @desc Renders block labels
 * @param {*} param0
 */
export const renderBlock = ({ model, className }) => `
<span class="block-icon-container">${model.get("icon")}</span>
<span>
${model.get("label")}
</span>
`
