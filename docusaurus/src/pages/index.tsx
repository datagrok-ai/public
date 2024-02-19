import React from 'react';
import clsx from 'clsx';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import Layout from '@theme/Layout';
import styles from './index.module.css';
import Hero from '../components/home/hero';
import Body from '../components/home/body';

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
      <HomepageHeader />
    </Layout>
  );
}
