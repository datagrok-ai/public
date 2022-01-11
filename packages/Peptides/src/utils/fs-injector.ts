export class FSInjector {
    stdinBuffer: string;
    stdoutBuffer: string;
    stderrBuffer: string;
    inIndex: number;

    constructor(stdinBuffer: string) {
        this.stdinBuffer = stdinBuffer;
        this.stdoutBuffer = '';
        this.stderrBuffer = '';
        this.inIndex = 0;
    }
    
    reset() {
        console.log('reset');
        this.stdoutBuffer = '';
        this.stderrBuffer = '';
        this.inIndex = 0;
    }

    stdin(): number | null {
        if (this.inIndex < this.stdinBuffer.length) {
          let code = this.stdinBuffer.charCodeAt(this.inIndex);
          ++this.inIndex;
          return code;
        } else {
          return null;
        }
      }
  
    stdout(code: number) {
        this.stdoutBuffer += String.fromCharCode(code);
        //console.log(this);
        //console.log('stdout');
    }
  
    stderr(code: number) {
        this.stderrBuffer += String.fromCharCode(code);
        //console.log(this.stderrBuffer);
    }
}