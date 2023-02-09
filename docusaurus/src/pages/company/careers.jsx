import Markdown from 'markdown-to-jsx';
import React, { useEffect, useState, Component } from 'react';

import Fonts from '/static/docusaurus_css/fonts.css';
import Icons from '/static/docusaurus_css/all.min.css';
import Bootstrap from '/static/docusaurus_css/bootstrap.min.css';
import SignupLogin from '/static/docusaurus_css/signup_login.css';
import Datagrok from '/static/docusaurus_css/datagrok.css';

import MainMenu from '@site/src/components/default/menu.jsx';
import CareersLayout from '@site/src/docs/careers/default.mdx';
import Banner from '@site/src/components/default/banner';
import Footer from '@site/src/components/default/footer';
  
export default function root() {
    return(
        <div className='bg-light'>
            <MainMenu name='Careers'/>
            <section className="container-fluid bg-light py-5">
                <div className="container">
                    <CareersLayout/>
                </div>
            </section>
            <Banner/>
            <Footer/>
        </div>
    )
}