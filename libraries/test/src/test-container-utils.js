var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import * as grok from 'datagrok-api/grok';
import { awaitCheck } from './test';
export function ensureContainerRunning(containerName, timeout = 300000) {
    return __awaiter(this, void 0, void 0, function* () {
        const container = yield grok.dapi.docker.dockerContainers.filter(containerName).first();
        if (!(container.status.startsWith('started') || container.status.startsWith('checking'))) {
            console.log(`starting container ${container.name}`);
            yield grok.dapi.docker.dockerContainers.run(container.id, false);
        }
        let started = false;
        yield awaitCheck(() => {
            grok.dapi.docker.dockerContainers.find(container.id).then((cont) => {
                started = cont.status.startsWith('started') || cont.status.startsWith('checking');
            });
            return started;
        }, `${containerName} hasn't been started after ${timeout / 60000} minutes`, timeout, 5000);
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC1jb250YWluZXItdXRpbHMuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyJ0ZXN0LWNvbnRhaW5lci11dGlscy50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7Ozs7QUFBQSxPQUFPLEtBQUssSUFBSSxNQUFNLG1CQUFtQixDQUFDO0FBQzFDLE9BQU8sRUFBRSxVQUFVLEVBQUUsTUFBTSxRQUFRLENBQUM7QUFFcEMsTUFBTSxVQUFnQixzQkFBc0IsQ0FBQyxhQUFxQixFQUFFLFVBQWtCLE1BQU07O1FBQzFGLE1BQU0sU0FBUyxHQUFHLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxNQUFNLENBQUMsZ0JBQWdCLENBQUMsTUFBTSxDQUFDLGFBQWEsQ0FBQyxDQUFDLEtBQUssRUFBRSxDQUFDO1FBQ3hGLElBQUksQ0FBQyxDQUFDLFNBQVMsQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFNBQVMsQ0FBQyxJQUFJLFNBQVMsQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFVBQVUsQ0FBQyxDQUFDLEVBQUU7WUFDeEYsT0FBTyxDQUFDLEdBQUcsQ0FBQyxzQkFBc0IsU0FBUyxDQUFDLElBQUksRUFBRSxDQUFDLENBQUM7WUFDcEQsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLE1BQU0sQ0FBQyxnQkFBZ0IsQ0FBQyxHQUFHLENBQUMsU0FBUyxDQUFDLEVBQUUsRUFBRSxLQUFLLENBQUMsQ0FBQztTQUNsRTtRQUVELElBQUksT0FBTyxHQUFHLEtBQUssQ0FBQztRQUNwQixNQUFNLFVBQVUsQ0FBQyxHQUFHLEVBQUU7WUFDcEIsSUFBSSxDQUFDLElBQUksQ0FBQyxNQUFNLENBQUMsZ0JBQWdCLENBQUMsSUFBSSxDQUFDLFNBQVMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxJQUFJLEVBQUUsRUFBRTtnQkFDakUsT0FBTyxHQUFHLElBQUksQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFNBQVMsQ0FBQyxJQUFJLElBQUksQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFVBQVUsQ0FBQyxDQUFDO1lBQ3BGLENBQUMsQ0FBQyxDQUFDO1lBQ0gsT0FBTyxPQUFPLENBQUM7UUFDakIsQ0FBQyxFQUFFLEdBQUcsYUFBYSw4QkFBOEIsT0FBTyxHQUFHLEtBQUssVUFBVSxFQUFFLE9BQU8sRUFBRSxJQUFJLENBQUMsQ0FBQztJQUM3RixDQUFDO0NBQUEiLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgKiBhcyBncm9rIGZyb20gJ2RhdGFncm9rLWFwaS9ncm9rJztcclxuaW1wb3J0IHsgYXdhaXRDaGVjayB9IGZyb20gJy4vdGVzdCc7XHJcblxyXG5leHBvcnQgYXN5bmMgZnVuY3Rpb24gZW5zdXJlQ29udGFpbmVyUnVubmluZyhjb250YWluZXJOYW1lOiBzdHJpbmcsIHRpbWVvdXQ6IG51bWJlciA9IDMwMDAwMCkge1xyXG4gIGNvbnN0IGNvbnRhaW5lciA9IGF3YWl0IGdyb2suZGFwaS5kb2NrZXIuZG9ja2VyQ29udGFpbmVycy5maWx0ZXIoY29udGFpbmVyTmFtZSkuZmlyc3QoKTtcclxuICBpZiAoIShjb250YWluZXIuc3RhdHVzLnN0YXJ0c1dpdGgoJ3N0YXJ0ZWQnKSB8fCBjb250YWluZXIuc3RhdHVzLnN0YXJ0c1dpdGgoJ2NoZWNraW5nJykpKSB7XHJcbiAgICBjb25zb2xlLmxvZyhgc3RhcnRpbmcgY29udGFpbmVyICR7Y29udGFpbmVyLm5hbWV9YCk7XHJcbiAgICBhd2FpdCBncm9rLmRhcGkuZG9ja2VyLmRvY2tlckNvbnRhaW5lcnMucnVuKGNvbnRhaW5lci5pZCwgZmFsc2UpO1xyXG4gIH1cclxuICBcclxuICBsZXQgc3RhcnRlZCA9IGZhbHNlO1xyXG4gIGF3YWl0IGF3YWl0Q2hlY2soKCkgPT4ge1xyXG4gICAgZ3Jvay5kYXBpLmRvY2tlci5kb2NrZXJDb250YWluZXJzLmZpbmQoY29udGFpbmVyLmlkKS50aGVuKChjb250KSA9PiB7XHJcbiAgICAgIHN0YXJ0ZWQgPSBjb250LnN0YXR1cy5zdGFydHNXaXRoKCdzdGFydGVkJykgfHwgY29udC5zdGF0dXMuc3RhcnRzV2l0aCgnY2hlY2tpbmcnKTtcclxuICAgIH0pO1xyXG4gICAgcmV0dXJuIHN0YXJ0ZWQ7XHJcbiAgfSwgYCR7Y29udGFpbmVyTmFtZX0gaGFzbid0IGJlZW4gc3RhcnRlZCBhZnRlciAke3RpbWVvdXQgLyA2MDAwMH0gbWludXRlc2AsIHRpbWVvdXQsIDUwMDApO1xyXG59Il19