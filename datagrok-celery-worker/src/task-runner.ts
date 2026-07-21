/** Executes one FuncCall end to end (mirrors DatagrokTask's task lifecycle in
 *  datagrok_task.py: input streaming -> run -> output streaming -> result CALL). */
import {FanoutPublisher} from './fanout-publisher';
import {FuncCall, FuncCallStatus} from './func-call';
import {PackageHost} from './package-host';
import {Const, PipeClient} from './pipe-client';
import {CancelledError, TaskContext} from './progress';
import {RunningTask} from './pidbox';
import {Settings} from './settings';
import {marshalInput, marshalOutput, validateCall} from './marshal';
import {logError, logInfo, logWarn} from './logger';

export class TaskRunner {
  private readonly settings: Settings;
  private readonly publisher: FanoutPublisher;
  private readonly host: PackageHost;
  private _current: RunningTask | null = null;

  constructor(settings: Settings, publisher: FanoutPublisher, host: PackageHost) {
    this.settings = settings;
    this.publisher = publisher;
    this.host = host;
  }

  get current(): RunningTask | null {
    return this._current;
  }

  async run(call: FuncCall): Promise<void> {
    logInfo(`Running ${call.funcName}`, call.id);
    const context = new TaskContext(call.id, this.publisher);
    this._current = {call: call, context: context};
    let pipe: PipeClient | null = null;
    call.status = FuncCallStatus.RUNNING;
    try {
      validateCall(call);
      const impl = await this.host.resolve(call);
      if (call.requiresPipe) {
        pipe = new PipeClient(`${this.settings.pipeUrl}/${call.id}`, {
          'x-member-name': `celery-${this.settings.celeryName}`,
          'authorization': this.settings.pipeKey,
        }, {
          messageTimeoutMs: this.settings.wsMessageTimeoutSeconds * 1000,
          paramTimeoutMs: this.settings.paramTimeoutMinutes * 60000,
          batchSize: call.binaryBatchSize,
        });
        await pipe.connect();
        context.pipe = pipe;
      }
      for (const param of call.inputParams)
        if (param.isStreamable)
          param.value = await pipe!.receiveParam(param.name);
      const dg = this.host.dg;
      for (const param of call.inputParams)
        marshalInput(param, dg);
      const args = call.inputParams.map((p) => p.value);
      context.install();
      let value: any;
      try {
        value = await impl(...args);
      }
      finally {
        context.uninstall();
      }
      if (context.cancelled)
        throw new CancelledError(call.id);
      const output = call.outputParams[0];
      if (output != null) {
        if (value == null)
          logWarn('Task returned null', call.id);
        const marshaled = marshalOutput(output, value, dg);
        if (output.isStreamable && marshaled.bytes != null)
          await pipe!.sendParam(marshaled.bytes, marshaled.tags!);
      }
      call.status = FuncCallStatus.COMPLETED;
      logInfo(`Completed ${call.funcName}`, call.id);
    }
    catch (e: any) {
      if (e instanceof CancelledError || context.cancelled) {
        call.status = FuncCallStatus.CANCELED;
        logInfo(`Canceled ${call.funcName}`, call.id);
      }
      else {
        call.status = FuncCallStatus.ERROR;
        call.errorMessage = e?.message ?? String(e);
        call.errorStackTrace = typeof e?.stack === 'string' ? e.stack : String(e);
        logError(`Task failed: ${call.errorMessage}`, call.id);
      }
    }
    finally {
      await this.publishResult(call, context, pipe);
      pipe?.close();
      this._current = null;
    }
  }

  /** Sends the result: 'CALL <json>' over the pipe for streamable calls, otherwise a
   *  fanout 'call' message. Datlas accepts BOTH channels for any call
   *  (call_queue_service.dart completes serverCalls from either), so a broken pipe
   *  falls back to AMQP instead of retrying the pipe like the python lib does. */
  private async publishResult(call: FuncCall, context: TaskContext, pipe: PipeClient | null): Promise<void> {
    if (context.cancelPublished)
      return; // pidbox already published the Canceled CALL
    try {
      if (pipe != null && pipe.isConnected) {
        pipe.sendText(`${Const.CALL} ${JSON.stringify(call.toJson())}`);
        return;
      }
    }
    catch (e: any) {
      logWarn(`CALL send over the pipe failed, falling back to AMQP: ${e?.message ?? e}`, call.id);
    }
    await this.publisher.publish(call.toJson(), call.id, 'call');
  }
}
