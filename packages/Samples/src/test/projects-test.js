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
import { category, test } from '@datagrok-libraries/utils/src/test';
category('Projects', () => {
    test('coffee_company', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('coffee_company');
    }));
    test('chem_demo', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('chem_demo');
    }), { skipReason: 'GROK-13612' });
    test('chem_tsne', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('chem_tsne');
    }), { skipReason: 'GROK-13621' });
    test('demo_notebooks', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('demo_notebooks');
    }));
    test('demog', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('demog');
    }));
    test('ecg', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('ecg');
    }));
    test('eeg', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('eeg');
    }));
    test('game-of-thrones', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('game_of_thrones');
    }));
    test('models', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('models');
    }));
    test('pa-income-by-county', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('pa_income_by_county');
    }));
    test('it-network', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('it_network');
    }));
    test('plates', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('plates');
    }));
    test('time-series-decomposition', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('time_series_decomposition');
    }));
    test('wells', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('wells');
    }));
    test('pic50', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('pic50');
    }));
    test('zbb', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('zbb');
    }));
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoicHJvamVjdHMtdGVzdC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzIjpbInByb2plY3RzLXRlc3QudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7Ozs7Ozs7O0FBQUEsT0FBTyxLQUFLLElBQUksTUFBTSxtQkFBbUIsQ0FBQztBQUMxQyxPQUFPLEVBQWdCLFFBQVEsRUFBaUIsSUFBSSxFQUFDLE1BQU0sb0NBQW9DLENBQUM7QUFHaEcsUUFBUSxDQUFDLFVBQVUsRUFBRSxHQUFHLEVBQUU7SUFDdEIsSUFBSSxDQUFDLGdCQUFnQixFQUFFLEdBQVMsRUFBRTtRQUM5QixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxnQkFBZ0IsQ0FBQyxDQUFBO0lBQ25ELENBQUMsQ0FBQSxDQUFDLENBQUM7SUFDSCxJQUFJLENBQUMsV0FBVyxFQUFFLEdBQVMsRUFBRTtRQUN6QixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxXQUFXLENBQUMsQ0FBQTtJQUM5QyxDQUFDLENBQUEsRUFBRSxFQUFDLFVBQVUsRUFBRSxZQUFZLEVBQUMsQ0FBQyxDQUFDO0lBQy9CLElBQUksQ0FBQyxXQUFXLEVBQUUsR0FBUyxFQUFFO1FBQ3pCLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLFdBQVcsQ0FBQyxDQUFBO0lBQzlDLENBQUMsQ0FBQSxFQUFFLEVBQUMsVUFBVSxFQUFFLFlBQVksRUFBQyxDQUFDLENBQUM7SUFDL0IsSUFBSSxDQUFDLGdCQUFnQixFQUFFLEdBQVMsRUFBRTtRQUM5QixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxnQkFBZ0IsQ0FBQyxDQUFBO0lBQ25ELENBQUMsQ0FBQSxDQUFDLENBQUM7SUFDSCxJQUFJLENBQUMsT0FBTyxFQUFFLEdBQVMsRUFBRTtRQUNyQixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsQ0FBQTtJQUMxQyxDQUFDLENBQUEsQ0FBQyxDQUFDO0lBQ0gsSUFBSSxDQUFDLEtBQUssRUFBRSxHQUFTLEVBQUU7UUFDbkIsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLFFBQVEsQ0FBQyxJQUFJLENBQUMsS0FBSyxDQUFDLENBQUE7SUFDeEMsQ0FBQyxDQUFBLENBQUMsQ0FBQztJQUNILElBQUksQ0FBQyxLQUFLLEVBQUUsR0FBUyxFQUFFO1FBQ25CLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLEtBQUssQ0FBQyxDQUFBO0lBQ3hDLENBQUMsQ0FBQSxDQUFDLENBQUM7SUFDSCxJQUFJLENBQUMsaUJBQWlCLEVBQUUsR0FBUyxFQUFFO1FBQy9CLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLGlCQUFpQixDQUFDLENBQUE7SUFDcEQsQ0FBQyxDQUFBLENBQUMsQ0FBQztJQUNILElBQUksQ0FBQyxRQUFRLEVBQUUsR0FBUyxFQUFFO1FBQ3RCLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLFFBQVEsQ0FBQyxDQUFBO0lBQzNDLENBQUMsQ0FBQSxDQUFDLENBQUM7SUFDSCxJQUFJLENBQUMscUJBQXFCLEVBQUUsR0FBUyxFQUFFO1FBQ25DLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLHFCQUFxQixDQUFDLENBQUE7SUFDeEQsQ0FBQyxDQUFBLENBQUMsQ0FBQztJQUNILElBQUksQ0FBQyxZQUFZLEVBQUUsR0FBUyxFQUFFO1FBQzFCLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLFlBQVksQ0FBQyxDQUFBO0lBQy9DLENBQUMsQ0FBQSxDQUFDLENBQUM7SUFDSCxJQUFJLENBQUMsUUFBUSxFQUFFLEdBQVMsRUFBRTtRQUN0QixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsQ0FBQTtJQUMzQyxDQUFDLENBQUEsQ0FBQyxDQUFDO0lBQ0gsSUFBSSxDQUFDLDJCQUEyQixFQUFFLEdBQVMsRUFBRTtRQUN6QyxNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQywyQkFBMkIsQ0FBQyxDQUFBO0lBQzlELENBQUMsQ0FBQSxDQUFDLENBQUM7SUFDSCxJQUFJLENBQUMsT0FBTyxFQUFFLEdBQVMsRUFBRTtRQUNyQixNQUFNLElBQUksQ0FBQyxJQUFJLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsQ0FBQTtJQUMxQyxDQUFDLENBQUEsQ0FBQyxDQUFDO0lBQ0gsSUFBSSxDQUFDLE9BQU8sRUFBRSxHQUFTLEVBQUU7UUFDckIsTUFBTSxJQUFJLENBQUMsSUFBSSxDQUFDLFFBQVEsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUE7SUFDMUMsQ0FBQyxDQUFBLENBQUMsQ0FBQztJQUNILElBQUksQ0FBQyxLQUFLLEVBQUUsR0FBUyxFQUFFO1FBQ25CLE1BQU0sSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLEtBQUssQ0FBQyxDQUFBO0lBQ3hDLENBQUMsQ0FBQSxDQUFDLENBQUM7QUFDUCxDQUFDLENBQUMsQ0FBQyIsInNvdXJjZXNDb250ZW50IjpbImltcG9ydCAqIGFzIGdyb2sgZnJvbSAnZGF0YWdyb2stYXBpL2dyb2snO1xuaW1wb3J0IHthZnRlciwgYmVmb3JlLCBjYXRlZ29yeSwgZGVsYXksIGV4cGVjdCwgdGVzdH0gZnJvbSAnQGRhdGFncm9rLWxpYnJhcmllcy91dGlscy9zcmMvdGVzdCc7XG5cblxuY2F0ZWdvcnkoJ1Byb2plY3RzJywgKCkgPT4ge1xuICAgIHRlc3QoJ2NvZmZlZV9jb21wYW55JywgYXN5bmMgKCkgPT4ge1xuICAgICAgICBhd2FpdCBncm9rLmRhcGkucHJvamVjdHMub3BlbignY29mZmVlX2NvbXBhbnknKVxuICAgIH0pO1xuICAgIHRlc3QoJ2NoZW1fZGVtbycsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ2NoZW1fZGVtbycpXG4gICAgfSwge3NraXBSZWFzb246ICdHUk9LLTEzNjEyJ30pO1xuICAgIHRlc3QoJ2NoZW1fdHNuZScsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ2NoZW1fdHNuZScpXG4gICAgfSwge3NraXBSZWFzb246ICdHUk9LLTEzNjIxJ30pO1xuICAgIHRlc3QoJ2RlbW9fbm90ZWJvb2tzJywgYXN5bmMgKCkgPT4ge1xuICAgICAgICBhd2FpdCBncm9rLmRhcGkucHJvamVjdHMub3BlbignZGVtb19ub3RlYm9va3MnKVxuICAgIH0pO1xuICAgIHRlc3QoJ2RlbW9nJywgYXN5bmMgKCkgPT4ge1xuICAgICAgICBhd2FpdCBncm9rLmRhcGkucHJvamVjdHMub3BlbignZGVtb2cnKVxuICAgIH0pO1xuICAgIHRlc3QoJ2VjZycsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ2VjZycpXG4gICAgfSk7XG4gICAgdGVzdCgnZWVnJywgYXN5bmMgKCkgPT4ge1xuICAgICAgICBhd2FpdCBncm9rLmRhcGkucHJvamVjdHMub3BlbignZWVnJylcbiAgICB9KTtcbiAgICB0ZXN0KCdnYW1lLW9mLXRocm9uZXMnLCBhc3luYyAoKSA9PiB7XG4gICAgICAgIGF3YWl0IGdyb2suZGFwaS5wcm9qZWN0cy5vcGVuKCdnYW1lX29mX3Rocm9uZXMnKVxuICAgIH0pO1xuICAgIHRlc3QoJ21vZGVscycsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ21vZGVscycpXG4gICAgfSk7XG4gICAgdGVzdCgncGEtaW5jb21lLWJ5LWNvdW50eScsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ3BhX2luY29tZV9ieV9jb3VudHknKVxuICAgIH0pO1xuICAgIHRlc3QoJ2l0LW5ldHdvcmsnLCBhc3luYyAoKSA9PiB7XG4gICAgICAgIGF3YWl0IGdyb2suZGFwaS5wcm9qZWN0cy5vcGVuKCdpdF9uZXR3b3JrJylcbiAgICB9KTtcbiAgICB0ZXN0KCdwbGF0ZXMnLCBhc3luYyAoKSA9PiB7XG4gICAgICAgIGF3YWl0IGdyb2suZGFwaS5wcm9qZWN0cy5vcGVuKCdwbGF0ZXMnKVxuICAgIH0pO1xuICAgIHRlc3QoJ3RpbWUtc2VyaWVzLWRlY29tcG9zaXRpb24nLCBhc3luYyAoKSA9PiB7XG4gICAgICAgIGF3YWl0IGdyb2suZGFwaS5wcm9qZWN0cy5vcGVuKCd0aW1lX3Nlcmllc19kZWNvbXBvc2l0aW9uJylcbiAgICB9KTtcbiAgICB0ZXN0KCd3ZWxscycsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ3dlbGxzJylcbiAgICB9KTtcbiAgICB0ZXN0KCdwaWM1MCcsIGFzeW5jICgpID0+IHtcbiAgICAgICAgYXdhaXQgZ3Jvay5kYXBpLnByb2plY3RzLm9wZW4oJ3BpYzUwJylcbiAgICB9KTtcbiAgICB0ZXN0KCd6YmInLCBhc3luYyAoKSA9PiB7XG4gICAgICAgIGF3YWl0IGdyb2suZGFwaS5wcm9qZWN0cy5vcGVuKCd6YmInKVxuICAgIH0pO1xufSk7XG4iXX0=