import React from 'react';
import styles from './styles.module.css';

import Fonts from '/static/docusaurus_css/fonts.css';
import Icons from '/static/docusaurus_css/all.min.css';
import Bootstrap from '/static/docusaurus_css/bootstrap.min.css';
import SignupLogin from '/static/docusaurus_css/signup_login.css';
import Datagrok from '/static/docusaurus_css/datagrok.css';

import MainMenu from '@site/src/components/new/navbar.jsx';
import Banner from '@site/src/components/default/banner';
import Footer from '@site/src/components/new/footer';

export default function VacancyPage({meta, children}) {
    return (
        <div className='bg-light'>
            <div className='new'><MainMenu name='Careers'/></div>
            <section className='container-fluid bg-light py-5'>
                <div className='container'>
                    <div className='row py-5'>
                        <div className='col-lg-8 col-md-10 col-sm-12 mx-auto'>
                            <a href='/company/careers' className='text-secondary d-inline-block mb-4'>← All positions</a>
                            <h1>{meta.name}</h1>
                            <div className='d-flex flex-row mb-4'>
                                <span className={'badge badge-pill mr-2 ' + styles.badgedg}>{meta.type}</span>
                                <span className={'badge badge-pill badge-info ' + styles.badgedg}>{meta.location}</span>
                            </div>
                            <div className='text-left'>
                                {children}
                            </div>
                            <br/>
                            <a className='btn btn-primary px-4' href={'mailto:hr@datagrok.ai?subject=' + meta.name}>Apply</a>
                        </div>
                    </div>
                </div>
            </section>
            <Banner/>
            <div className='new'><Footer/></div>
        </div>
    );
}
