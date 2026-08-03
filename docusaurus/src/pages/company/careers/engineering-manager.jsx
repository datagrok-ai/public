import React from 'react';
import VacancyPage from '@site/src/components/careers/vacancy-page.jsx';
import Content, {meta} from '@site/src/docs/careers/engineering-manager.mdx';

export default function Page() {
    return (
        <VacancyPage meta={meta}>
            <Content/>
        </VacancyPage>
    );
}
