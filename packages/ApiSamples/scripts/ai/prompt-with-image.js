// An image attachment: a Blob whose type is the MIME type. Engines that cannot take images throw.

const bytes = await grok.dapi.files.readAsBytes('System:DemoFiles/images/cats/dog1.jpg');
const image = new Blob([bytes], {type: 'image/jpeg'});
const {text} = await grok.ai.prompt('What animal is this?', {attachments: [{type: 'image', data: image}]});
grok.shell.info(text);
