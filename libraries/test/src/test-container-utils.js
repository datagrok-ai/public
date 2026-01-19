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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoidGVzdC1jb250YWluZXItdXRpbHMuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyJ0ZXN0LWNvbnRhaW5lci11dGlscy50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7Ozs7Ozs7QUFBQSxPQUFRLEtBQUssSUFBSSxNQUFNLG1CQUFtQixDQUFDO0FBQzNDLE9BQU8sRUFBRSxVQUFVLEVBQUUsTUFBTSxRQUFRLENBQUM7QUFFcEMsTUFBTSxVQUFnQixzQkFBc0IsQ0FBQyxhQUFxQixFQUFFLFVBQWtCLE1BQU07O1FBQzFGLE1BQU0sU0FBUyxHQUFHLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxNQUFNLENBQUMsZ0JBQWdCLENBQUMsTUFBTSxDQUFDLGFBQWEsQ0FBQyxDQUFDLEtBQUssRUFBRSxDQUFDO1FBQ3hGLElBQUksQ0FBQyxDQUFDLFNBQVMsQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFNBQVMsQ0FBQyxJQUFJLFNBQVMsQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFVBQVUsQ0FBQyxDQUFDLEVBQUU7WUFDeEYsT0FBTyxDQUFDLEdBQUcsQ0FBQyxzQkFBc0IsU0FBUyxDQUFDLElBQUksRUFBRSxDQUFDLENBQUM7WUFDcEQsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLE1BQU0sQ0FBQyxnQkFBZ0IsQ0FBQyxHQUFHLENBQUMsU0FBUyxDQUFDLEVBQUUsRUFBRSxLQUFLLENBQUMsQ0FBQztTQUNsRTtRQUVELElBQUksT0FBTyxHQUFHLEtBQUssQ0FBQztRQUNwQixNQUFNLFVBQVUsQ0FBQyxHQUFHLEVBQUU7WUFDcEIsSUFBSSxDQUFDLElBQUksQ0FBQyxNQUFNLENBQUMsZ0JBQWdCLENBQUMsSUFBSSxDQUFDLFNBQVMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxJQUFJLEVBQUUsRUFBRTtnQkFDakUsT0FBTyxHQUFHLElBQUksQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFNBQVMsQ0FBQyxJQUFJLElBQUksQ0FBQyxNQUFNLENBQUMsVUFBVSxDQUFDLFVBQVUsQ0FBQyxDQUFDO1lBQ3BGLENBQUMsQ0FBQyxDQUFDO1lBQ0gsT0FBTyxPQUFPLENBQUM7UUFDakIsQ0FBQyxFQUFFLEdBQUcsYUFBYSw4QkFBOEIsT0FBTyxHQUFHLEtBQUssVUFBVSxFQUFFLE9BQU8sRUFBRSxJQUFJLENBQUMsQ0FBQztJQUM3RixDQUFDO0NBQUEiLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgICogYXMgZ3JvayBmcm9tICdkYXRhZ3Jvay1hcGkvZ3Jvayc7XHJcbmltcG9ydCB7IGF3YWl0Q2hlY2sgfSBmcm9tICcuL3Rlc3QnO1xyXG5cclxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIGVuc3VyZUNvbnRhaW5lclJ1bm5pbmcoY29udGFpbmVyTmFtZTogc3RyaW5nLCB0aW1lb3V0OiBudW1iZXIgPSAzMDAwMDApIHtcclxuICBjb25zdCBjb250YWluZXIgPSBhd2FpdCBncm9rLmRhcGkuZG9ja2VyLmRvY2tlckNvbnRhaW5lcnMuZmlsdGVyKGNvbnRhaW5lck5hbWUpLmZpcnN0KCk7XHJcbiAgaWYgKCEoY29udGFpbmVyLnN0YXR1cy5zdGFydHNXaXRoKCdzdGFydGVkJykgfHwgY29udGFpbmVyLnN0YXR1cy5zdGFydHNXaXRoKCdjaGVja2luZycpKSkge1xyXG4gICAgY29uc29sZS5sb2coYHN0YXJ0aW5nIGNvbnRhaW5lciAke2NvbnRhaW5lci5uYW1lfWApO1xyXG4gICAgYXdhaXQgZ3Jvay5kYXBpLmRvY2tlci5kb2NrZXJDb250YWluZXJzLnJ1bihjb250YWluZXIuaWQsIGZhbHNlKTtcclxuICB9XHJcbiAgXHJcbiAgbGV0IHN0YXJ0ZWQgPSBmYWxzZTtcclxuICBhd2FpdCBhd2FpdENoZWNrKCgpID0+IHtcclxuICAgIGdyb2suZGFwaS5kb2NrZXIuZG9ja2VyQ29udGFpbmVycy5maW5kKGNvbnRhaW5lci5pZCkudGhlbigoY29udCkgPT4ge1xyXG4gICAgICBzdGFydGVkID0gY29udC5zdGF0dXMuc3RhcnRzV2l0aCgnc3RhcnRlZCcpIHx8IGNvbnQuc3RhdHVzLnN0YXJ0c1dpdGgoJ2NoZWNraW5nJyk7XHJcbiAgICB9KTtcclxuICAgIHJldHVybiBzdGFydGVkO1xyXG4gIH0sIGAke2NvbnRhaW5lck5hbWV9IGhhc24ndCBiZWVuIHN0YXJ0ZWQgYWZ0ZXIgJHt0aW1lb3V0IC8gNjAwMDB9IG1pbnV0ZXNgLCB0aW1lb3V0LCA1MDAwKTtcclxufSJdfQ==