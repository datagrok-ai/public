import React from "react"
import JSONData from "./plugins.json"

function makeid(length) {
    let result = '';
    const characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    const charactersLength = characters.length;
    let counter = 0;
    while (counter < length) {
        result += characters.charAt(Math.floor(Math.random() * charactersLength));
        counter += 1;
    }
    return result;
}

function textToLink(text) {
    let parts = text.split(/(\[[\S\s]+\]\(\S+\))/);
    for (let i = 1; i < parts.length; i += 2) {
        let regex = /\[([\S\s]+)\]\((\S+)\)/;
        const groups = regex.exec(parts[i]);
        const id = makeid(5);
        if (regex.test(parts[i])) {
            const groups = regex.exec(parts[i]);
            parts[i] = <a key={'link' + id}
                          href={groups[2].replace(/(\.\.\/)+help/, '/help').replace(/\.md/, '')}>{groups[1]}</a>;
        } else {
            parts[i] = parts[i] + ' ';
        }
    }
    return parts;
}

const JSONTable = () => (<table>
    <thead>
    <tr>
        {['Name', 'Date', 'Version', 'Summary'].map((header, index) => (<th key={'header' + index} style={{
            width: index === 1 ? 112 : "auto",
        }}>{header}</th>))}
    </tr>
    </thead>
    <tbody>
    {JSONData.map((row, rowIndex) => (<tr key={'row' + rowIndex}>
        <td key='cell0'>{row.readme ? <a key={'link' + makeid(5)}
                                         href={"https://github.com/datagrok-ai/public/blob/master/" + row.dir + "/README.md"}>{row.name}</a> : row.name}
        </td>
        <td key='cell1'>{row.changelogDate}</td>
        <td key='cell2'>{row.changelog ? <a key={'link' + makeid(5)}
                                            href={"/help/deploy/releases/plugins/" + row.name.replaceAll(' ', '_') + "#" + row.version.replaceAll('.', '') + '-' + row.changelogDate}>{row.version}</a> : row.version}</td>
        <td key='cell3'>{row.summary.split('\n').map((line) => <div key={'div' + makeid(5)}>
            {textToLink(line.replaceAll(/### ([A-z ]+)/g, '$1:'))}
        </div>)}</td>
    </tr>))}
    </tbody>
</table>)

// JSONData.map((row, rowIndex) => row.summary.split('\n').map((line) => console.log(textToLink(line.replaceAll(/### ([A-z ]+)/g, '$1:')))))

export default JSONTable
