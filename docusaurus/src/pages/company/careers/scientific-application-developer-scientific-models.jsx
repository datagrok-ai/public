import React from 'react';
import VacancyPage from '@site/src/components/careers/vacancy-page.jsx';
import Content, {meta} from '@site/src/docs/careers/scientific-application-developer-model-hub.mdx';

export default function Page() {
    return (
        <VacancyPage meta={meta}>
            <Content/>
        </VacancyPage>
    );
}
