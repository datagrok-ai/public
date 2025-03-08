import React from 'react';
import clsx from 'clsx';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import Layout from '@theme/Layout';
import styles from './index.module.css';
import Hero from '../components/home/hero';
import Body from '../components/home/body';
import Head from '@docusaurus/Head';

function HomepageHeader() {
  const {siteConfig} = useDocusaurusContext();
  // window.location.href = 'https://datagrok.ai';
  return (
    <div>
      <header>
      <Hero/>
      <Body/>
    </header>
    </div>
  );
}

export default function Home(): JSX.Element {
  const {siteConfig} = useDocusaurusContext();
  return (
    <Layout
      title={`${siteConfig.title}`}>
        <Head>
            <script async src="https://www.googletagmanager.com/gtag/js?id=G-KGZGSFPGHS"></script>
            <script>
                window.dataLayer = window.dataLayer || [];
                function gtag(){window.dataLayer.push(arguments)}
                gtag('js', new Date());

                gtag('config', 'G-KGZGSFPGHS');
            </script>
        </Head>
      <HomepageHeader />
    </Layout>
  );
}
