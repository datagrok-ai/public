import * as grok from 'datagrok-api/grok';

export function removeChildren(node) {
    while (node.firstChild)
        node.removeChild(node.firstChild);
}

export async function setupEnvironment(environment, containerId) {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, `/notebook/helper/environments/${environment.name}/exists`,
        {method: 'POST', body: environment.environment});
    if (response.status !== 200 && response.status !== 404) {
        const body = await response.text()
        throw body ?? response.statusText;
    }
    if (response.status === 200)
        return;
    const createResponse = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, `/notebook/helper/environments/${environment.name}/create`,
        {method: 'POST', body: environment.environment});
    const createBody = await createResponse.text();
    if (createBody !== 'created')
        throw createBody;
}

export async function editNotebook(notebook, containerId) {
    const resp = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, '/notebook/api');
    if (resp.status !== 200 && resp.status !== 201) {
        const error = (await resp.json()).error;
        throw new Error(`Jupyter Notebook is not available: ${error}`);
    }
    const data = notebook.notebook;
    if (!data['metadata']['datagrok'])
        data['metadata']['datagrok'] = {};
    const metadata = data['metadata']['datagrok'];
    metadata["id"] = notebook.id;
    metadata["host"] = getApiRoot();
    metadata["user_id"] = (await grok.dapi.users.current()).id;
    metadata["session_token"] = getAuthToken();
    notebook.notebook = data;
    const fileName = `${crypto.randomUUID()}.ipynb`;
    return await createFile(containerId,'', fileName, JSON.stringify(data));
}

async function createFile(containerId, path, name, content = "") {
    const filePath = (path === '') ? name : `${path}/${name}`;
    const params = (content === '') ? {"type": "file"} : {"type": "file", "format": "text", "content": content};
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, `/notebook/api/contents/${filePath}`,
        {method: 'PUT', body: JSON.stringify(params), headers: {"Content-Type": "application/json"}});
    if (response.status > 201)
        throw response.statusText;
    return filePath;
}

function getApiRoot() {
    let parsedUrl = new URL(window.location.href);
    return `${parsedUrl.protocol}//${parsedUrl.host}/api`;
}

export function getAuthToken() {
    let cookies = document.cookie.split("; ");
    for (let cookie of cookies) {
        let [key, value] = cookie.split("=");
        if (key === 'auth') return decodeURIComponent(value);
    }
    return null;
}
