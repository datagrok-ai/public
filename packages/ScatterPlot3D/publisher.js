const http = require('http');
//const { exec } = require('node:child_process');
const webs = require('ws')
const { exec } = require('child_process');

var ws = new webs('http://10.0.75.1:8024');

ws.on('open', () => {
    console.log('conn');
    ws.send('heelo from container')
})

ws.onmessage = function (event) {
    console.log('event.data: ',event.data);
    exec('webpack-cli'
    , (err, stdout, stderr) => {
        if (err) {
          // node couldn't execute the command
          console.log('if err: ', err)
          return;
        }
      
        // the *entire* stdout and stderr (buffered)
        console.log(Date.now())
        console.log(`stdout: ${stdout}`);
        console.log(`stderr: ${stderr}`);
        ws.send('webpacked1')
        exec('grok publish https://dev.datagrok.ai/api --key 952e8090-7663-53c0-a49a-81ade66950be --rebuild'
          , (err2, stdout2, stderr2) => {
            if (err2) {
              console.error('if err2 ', err2)
              return 1;
            }
                    // the *entire* stdout and stderr (buffered)
          console.log(Date.now())
          console.log(`stdout: ${stdout}`);
          console.log(`stderr: ${stderr}`);
        ws.send('compiled1')

          } // (err2, stdout2, strerr2)
        )

      })
    
    

  }