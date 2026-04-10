import React from 'react';
import VacancyPage from '@site/src/components/careers/vacancy-page.jsx';
import Content, {meta} from '@site/src/docs/careers/software-dev-dd.mdx';

export default function Page() {
    return (
        <VacancyPage meta={meta}>
            <Content/>
        </VacancyPage>
    );
}
