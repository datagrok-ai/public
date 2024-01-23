import React from 'react';
import Fonts from '/static/docusaurus_css/fonts.css';
import Icons from '/static/docusaurus_css/all.min.css';
import Animate from '/static/docusaurus_css/new/animate.min.css';
import MainCSS from '/static/docusaurus_css/new/main.css';
import Head from '@docusaurus/Head';
import Layout from '@theme/Layout/Provider';
import NavBar from '../../components/new/navbar';
import Footer from '../../components/new/footer';
import HomeLayout from '../../docs/home.mdx';

export default function root() {
    return(
        <>
        <Layout>
        <Head>
            <meta name="robots" content="noindex" />
            <meta name="googlebot" content="noindex" />
            <title>Datagrok</title>
            <meta name="description" content="Built for any data, tailor-made for life sciences" />
        </Head>
        <div className='new'>
            <NavBar page={''}/>
            <div class="container-fluid text-center px-0">
            <HomeLayout/>
            </div>
            <Footer/>
        </div> 
        </Layout>
        </>
    );
}


