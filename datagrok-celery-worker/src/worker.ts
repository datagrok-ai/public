#!/usr/bin/env node
/** Entrypoint: wires settings, health endpoint, pidbox + task queue consumers, and
 *  graceful shutdown. */
import {CeleryConsumer} from './celery-consumer';
import {FanoutPublisher} from './fanout-publisher';
import {PackageHost} from './package-host';
import {PidboxConsumer, RevokedSet} from './pidbox';
import {Settings} from './settings';
import {TaskRunner} from './task-runner';
import {startHealthServer} from './health';
import {logError, logInfo} from './logger';

const SHUTDOWN_WAIT_MS = 30000;

async function main(): Promise<void> {
  let settings: Settings;
  try {
    settings = Settings.fromEnv();
  }
  catch (e: any) {
    console.error(`datagrok-celery-worker: ${e?.message ?? e}`);
    process.exit(1);
    return;
  }

  logInfo(`Starting datagrok-celery-worker for package '${settings.packageName}'` +
    `${settings.packageVersion ? ` v${settings.packageVersion}` : ''}` +
    `, queue '${settings.taskQueueName}', hostname '${settings.celeryHostname}'`);

  const health = startHealthServer(settings.healthPort);
  const publisher = new FanoutPublisher(settings.brokerUrl);
  const revoked = new RevokedSet();
  const host = new PackageHost(settings);
  const runner = new TaskRunner(settings, publisher, host);
  const pidbox = new PidboxConsumer(settings, publisher, revoked, () => runner.current);
  const consumer = new CeleryConsumer(settings, publisher, revoked, (call) => runner.run(call));

  await pidbox.start();
  await consumer.start();

  let shuttingDown = false;
  const shutdown = async (signal: string): Promise<void> => {
    if (shuttingDown)
      return;
    shuttingDown = true;
    logInfo(`Received ${signal}, shutting down...`);
    try {
      await consumer.stopConsuming();
      if (!await consumer.waitForIdle(SHUTDOWN_WAIT_MS))
        logError(`In-flight task did not finish within ${SHUTDOWN_WAIT_MS} ms, exiting anyway`);
      await consumer.close();
      await pidbox.close();
      await publisher.close();
      health.close();
    }
    catch (e: any) {
      logError(`Shutdown error: ${e?.message ?? e}`);
    }
    process.exit(0);
  };
  process.on('SIGTERM', () => void shutdown('SIGTERM'));
  process.on('SIGINT', () => void shutdown('SIGINT'));
}

main().catch((e) => {
  logError(`Fatal: ${e?.message ?? e}`);
  process.exit(1);
});
