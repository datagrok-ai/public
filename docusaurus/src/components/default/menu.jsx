import styles from './styles.module.css';
import React from 'react';

export default function mainMenu (page){
    console.log(page.name);
    return(
        <nav className="navbar fixed-top navbar-expand-lg navbar-light bg-white">
            <a className="navbar-brand mb-0 h1" href="https://datagrok.ai/">datagrok</a>
            <button className="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbar" aria-controls=""
                aria-expanded="false" aria-label="Toggle navigation">
                <span className="navbar-toggler-icon"></span>
            </button>

            <div className="collapse w-100 navbar-collapse" id="navbar">
                <ul className="navbar-nav ml-auto mb-4 my-lg-0 text-uppercase">
                    <li className="nav-item dropdown">
                        <a className="nav-link dropdown-toggle" role="button" id="navbarFeatures" data-toggle="dropdown"
                            href="https://datagrok.ai/#">Features</a>
                        <div className="dropdown-menu" aria-labelledby="navbarFeatures">
                            <a className="dropdown-item" href="https://datagrok.ai/#access">Access</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#govern">Govern</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#transform">Transform</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#explore">Explore</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#learn">Learn</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#share">Share</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#apply">Apply</a>
                            <a className="dropdown-item" href="https://datagrok.ai/#deploy">Deploy &amp; Integrate</a>
                        </div>
                    </li>

                    <li className="nav-item dropdown">
                        <a className="nav-link dropdown-toggle" role="button" id="navbarSolutions" data-toggle="dropdown"
                            href="https://datagrok.ai/#" aria-expanded="false">Solutions</a>
                        <div className="dropdown-menu" aria-labelledby="navbarSolutions" id="solutionsMenu">
                            <a className="dropdown-item" href="https://datagrok.ai/">Datagrok</a>
                            <a className="dropdown-item" href="https://datagrok.ai/enterprise">Enterprise</a>
                            <a className="dropdown-item" href="https://datagrok.ai/data-science">Data Science</a>
                            <a className="dropdown-item" href="https://datagrok.ai/cheminformatics">Cheminformatics</a>
                            <div className="dropdown-divider"></div>
                            <a className="dropdown-item" href="https://datagrok.ai/solutions">Customer Use Cases</a>
                        </div>
                    </li>
                    <li className="nav-item dropdown">
                        <a className="nav-link dropdown-toggle" role="button" id="navbarCompany" data-toggle="dropdown"
                            href="https://datagrok.ai/#" aria-expanded="false">Company</a>
                        <div className="dropdown-menu" aria-labelledby="navbarCompany">
                        <a className={page.name == 'Careers' ? "dropdown-item active":"dropdown-item"} href="https://datagrok.ai/company/careers">Careers</a>
                            <a className={page.name == 'Team' ? "dropdown-item active":"dropdown-item"} href="https://datagrok.ai/company/team">Team</a>
                        </div>
                    </li>
                    <li className="nav-item dropdown">
                        <a className="nav-link dropdown-toggle" role="button" id="navbarHelp" data-toggle="dropdown"
                            href="https://datagrok.ai/#">Help</a>
                        <div className="dropdown-menu" aria-labelledby="navbarHelp">
                            <a className="dropdown-item" href="https://datagrok.ai/help">Wiki</a>
                            <a className="dropdown-item" href="https://datagrok.ai/js-api">JS API</a>
                            <a className="dropdown-item" href="https://community.datagrok.ai/">Community</a>
                            <a className="dropdown-item" href="https://us02web.zoom.us/meeting/register/up0vfuCrpjgpHNzi371YEJIQ4GkMpTm4disW#/registration">User group meetings</a>
                            <a className="dropdown-item d4-link-external"
                                href="https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg"><i className="fab fa-youtube"
                                ></i>Youtube</a>
                        </div>
                    </li>
                    <li className="nav-item">
                        <a className="btn btn-primary px-4" target="_blank" href="https://public.datagrok.ai/">Launch</a>
                    </li>
                </ul>
            </div>
        </nav>
    );
}