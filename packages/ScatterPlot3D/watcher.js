const chok = require('chokidar');
const http = require('http');
const webs = require('ws')
//const { pup } = require('puppetier')
//const pup = require('puppeteer');


var wss = new webs.Server({
    port: 8024,

})

var a = {}
wss.on('connection', (ws) => {
    a.ws = ws;
    ws.on('message', (m) => {
        console.log('mm', m)
        if (m == 'compiled') {

        }
    } 
)})

chok.watch(['sf3/*.js', 'sf3/webgl/*.js', 'src/*.js', '*.js']).on('change', async (e,p) => {
	console.log('change event2: ', e);
    a.ws.send('saved2 !')
console.log(p);
return 0
let opts = { hostname: '172.17.248.243', port: 3000, method: 'GET'};
let req = http.request(opts, (res) => {
	console.log('responce');
});
req.end();
}); // watch --disable-gpu
