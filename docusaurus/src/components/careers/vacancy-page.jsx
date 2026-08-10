import React, {useState} from 'react';
import styles from './styles.module.css';

import Fonts from '/static/docusaurus_css/fonts.css';
import Icons from '/static/docusaurus_css/all.min.css';
import Bootstrap from '/static/docusaurus_css/bootstrap.min.css';
import SignupLogin from '/static/docusaurus_css/signup_login.css';
import Datagrok from '/static/docusaurus_css/datagrok.css';

import MainMenu from '@site/src/components/new/navbar.jsx';
import Banner from '@site/src/components/default/banner';
import Footer from '@site/src/components/new/footer';

const hrEmail = 'hr@datagrok.ai';

export default function VacancyPage({meta, children}) {
    const [copied, setCopied] = useState(false);

    const copyEmail = () => {
        navigator.clipboard.writeText(hrEmail);
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
    };

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
                            <a className='btn btn-primary px-4' href={'mailto:' + hrEmail + '?subject=' + encodeURIComponent(meta.name)}
                               data-toggle='modal' data-target='#applyDialog'>Apply</a>
                        </div>
                    </div>
                </div>
            </section>
            <div className='modal fade' tabIndex='-1' id='applyDialog' aria-hidden='true'>
                <div className='modal-dialog modal-dialog-centered'>
                    <div className='modal-content'>
                        <div className={'modal-body p-4 ' + styles.dialogBody}>
                            <button type='button' className='close' data-dismiss='modal' aria-label='Close'>
                                <span aria-hidden='true'>&times;</span>
                            </button>
                            <h4 className='mt-2'>We'd love to hear from you</h4>
                            <p className='mb-0'>
                                To apply for the <b>{meta.name}</b> position, drop us a line at{' '}
                                <span className='text-nowrap'>
                                    <a href={'mailto:' + hrEmail + '?subject=' + encodeURIComponent(meta.name)}>{hrEmail}</a>
                                    <i className={copied ? 'fas fa-check ' + styles.copyIconDone : 'far fa-copy ' + styles.copyIcon}
                                       title='Copy email to clipboard' onClick={copyEmail}/>
                                </span>
                                {' '}— tell us a bit about yourself and attach your CV.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
            <Banner/>
            <div className='new'><Footer/></div>
        </div>
    );
}
