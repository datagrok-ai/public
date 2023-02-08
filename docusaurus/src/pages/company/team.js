import React from 'react';
import Fonts from '/static/docusaurus_css/fonts.css';
import Icons from '/static/docusaurus_css/all.min.css';
import Bootstrap from '/static/docusaurus_css/bootstrap.min.css';
import SignupLogin from '/static/docusaurus_css/signup_login.css';
import Datagrok from '/static/docusaurus_css/datagrok.css';

function menu() {
    return (
        <nav class="navbar fixed-top navbar-expand-lg navbar-light bg-white">
            <a class="navbar-brand mb-0 h1" href="https://datagrok.ai/">datagrok</a>
            <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbar" aria-controls=""
                aria-expanded="false" aria-label="Toggle navigation">
                <span class="navbar-toggler-icon"></span>
            </button>

            <div class="collapse w-100 navbar-collapse" id="navbar">
                <ul class="navbar-nav ml-auto mb-4 my-lg-0 text-uppercase">
                    <li class="nav-item dropdown">
                        <a class="nav-link dropdown-toggle" role="button" id="navbarFeatures" data-toggle="dropdown"
                            href="https://datagrok.ai/#">Features</a>
                        <div class="dropdown-menu" aria-labelledby="navbarFeatures">
                            <a class="dropdown-item" href="https://datagrok.ai/#access">Access</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#govern">Govern</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#transform">Transform</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#explore">Explore</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#learn">Learn</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#share">Share</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#apply">Apply</a>
                            <a class="dropdown-item" href="https://datagrok.ai/#deploy">Deploy &amp; Integrate</a>
                        </div>
                    </li>

                    <li class="nav-item dropdown">
                        <a class="nav-link dropdown-toggle" role="button" id="navbarSolutions" data-toggle="dropdown"
                            href="https://datagrok.ai/#" aria-expanded="false">Solutions</a>
                        <div class="dropdown-menu" aria-labelledby="navbarSolutions" id="solutionsMenu">
                            <a class="dropdown-item" href="https://datagrok.ai/">Datagrok</a>
                            <a class="dropdown-item" href="https://datagrok.ai/enterprise">Enterprise</a>
                            <a class="dropdown-item" href="https://datagrok.ai/data-science">Data Science</a>
                            <a class="dropdown-item" href="https://datagrok.ai/cheminformatics">Cheminformatics</a>
                            <div class="dropdown-divider"></div>
                            <a class="dropdown-item" href="https://datagrok.ai/solutions">Customer Use Cases</a>
                        </div>
                    </li>
                    <li class="nav-item dropdown">
                        <a class="nav-link dropdown-toggle" role="button" id="navbarCompany" data-toggle="dropdown"
                            href="https://datagrok.ai/#" aria-expanded="false">Company</a>
                        <div class="dropdown-menu" aria-labelledby="navbarCompany">
                            <a class="dropdown-item" href="https://datagrok.ai/company/careers">Careers</a>
                            <a class="dropdown-item active" href="https://datagrok.ai/company/team">Team</a>
                        </div>
                    </li>
                    <li class="nav-item dropdown">
                        <a class="nav-link dropdown-toggle" role="button" id="navbarHelp" data-toggle="dropdown"
                            href="https://datagrok.ai/#">Help</a>
                        <div class="dropdown-menu" aria-labelledby="navbarHelp">
                            <a class="dropdown-item" href="https://datagrok.ai/help">Wiki</a>
                            <a class="dropdown-item" href="https://datagrok.ai/js-api">JS API</a>
                            <a class="dropdown-item" href="https://community.datagrok.ai/">Community</a>
                            <a class="dropdown-item d4-link-external"
                                href="https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg"><i class="fab fa-youtube"
                                ></i>Youtube</a>
                        </div>
                    </li>
                    <li class="nav-item">
                        <a class="btn btn-primary px-4" target="_blank" href="https://public.datagrok.ai/">Launch</a>
                    </li>
                </ul>
            </div>
        </nav>
    );
}

function hero() {
    return (
        <div class="row py-5 my-3 align-items-center">
            <div class="col-lg-7 mb-4 px-5">
                <p class="display-5 text-secondary">We work together to build great things and slove real-world problems that impact</p>
                <a href="https://datagrok.ai/company/careers" class="mt-3 btn btn-primary px-">JOIN OUR TEAM</a>
            </div>
            <div class="col-lg-5 my-3">
                <video width="100%" height="280" class="rounded pt-2" controls>
                    Your browser does not support the video tag.
                </video>
            </div>
        </div>
    );
}

function teamHero() {
    return (
        <div class="row py-3 align-items-center">
            <div class="col-lg-12">
                <h1 class="section-white mb-2">Meet the team</h1>
                <p class="text-secondary">A remote-first team of experts with the HQ in Philadelphia, PA.</p>
            </div>
        </div>
    );
}

function rednderCard(props) {
    return (
        <div class="col-lg-3 col-md-3 col-sm-6 col-6 my-3">
                <div class="card hover-shadow border-0 border-0 border-light shadow-sm h-100">
                    <a class="btn link" data-toggle="modal" data-target={"#"+props.id}>
                        <img src={props.avatar} class="rounded-circle mt-4" height="96" width="96"></img>
                        <div class="card-body pt-0">
                            <p class="my-3">
                                <p class="h6 text-dark mt-3 mb-1">{props.name}</p>
                                <p class="font-size-sm mb-2" style={{ color: 'var(--grey-4)' }}>{props.position}</p>
                            </p>
                        </div>
                    </a>
                </div>

            <div class="modal fade" tabindex="-1" id={props.id} aria-hidden="true">
                <div class="modal-dialog modal-dialog-scrollable">
                    <div class="modal-content">
                        <div class="modal-body text-center ">
                            <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                                <i class="far fa-times fa-xs text-secondary"></i>
                            </button>
                            <img src={props.avatar} class="rounded-circle mt-4 " height="124" width="124"></img>
                                <h4 class="modal-title mt-2" id="exampleModalLabel">{props.name}</h4>
                                <p class="font-size-sm" style={{color:'var(--blue-2)'}}>{props.position}</p>
                                {props.linkedin != undefined ? <a href={"https://www.linkedin.com/in/"+props.linkedin} target="_blank"><i class="fab fa-linkedin-in fa-lg mx-1"></i></a>:null}
                                {props.github != undefined ? <a href={"https://github.com/"+props.github} target="_blank"><i class="fab fa-github fa-lg mx-1"></i></a>:null}
                                <div class="text-left p-3">
                                    <div class="text-secondary mb-3 text-justify">
                                     {props.description}
                                    </div>
                                </div>

                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}

function contactus(){
    return (
        <div class="contact-us-wrapper">
            <div class="contact-us-button contact-us">request a demo</div>
        </div>
    );
}

function footer(){
    return (
        <footer class="py-3 px-4 text-light footer-landing">
            <div class="container-fluid">
                <div class="row">
                    <div class="col-md-6"></div>
                    <div class="col-md-6 text-right">
                        Datagrok, Inc<br/>
                            <a class="text-light" href="mailto:info@datagrok.ai">info@datagrok.ai</a>
                    </div>
                </div>
            </div>
        </footer>
    );
}

export default function root() {
    let props = {};
    return (
        <div class="root">
            {menu()}
                <section class="container-fluid bg-light py-5">
                    <div class="container">
                        <div class="row">
                            {hero()}
                        </div>
                    </div>
                </section>
                <section class="container-fluid bg-white my-5">
                    <div class="container bg-map">
                        <div class="row justify-content-md-center text-center">
                            {teamHero()}
                        </div>
                        <div class="row">
                            {rednderCard(props = {
                                name: 'Andrew Skalkin',
                                id: 'andrew-skalkin',
                                position: 'Chief Executive Officer',
                                avatar: '/docusaurus_img/team/as.jpeg',
                                description: `Software architect, data scientist, and entrepreneur with broad experience in building analytical solutions that apply techniques of data management, high-performance computing, data analysis, modeling, machine learning, usability, and efficient data visualization to enable knowledge discovery and decision support. Proven ability to deliver elegant and simple solutions to complex problems.
                                Past work includes architecting platforms in the areas of exploratory data analysis, data integration, data management, chemoinformatics, health informatics, LIMS, bioinformatics, biomarker research, drug development and manufacturing, statistical process control, and biosensor data integration.
                                Passionate about data science, performance, usability, visualizations, algorithms.`,
                                linkedin: 'andrew-skalkin',
                                github:'skalkin'
                            })}

                            {rednderCard(props = {
                                name: 'Larisa Bankurova',
                                id: 'larisa-bankurova',
                                position: 'Chief Financial Officer',
                                avatar: '/docusaurus_img/team/lb.jpeg',
                                description: 'Specializing in working with high-growth technology companies operating in the global marketplace. Strong competency in US GAAP and IFRS; hands-on IPO experience; private, pre-IPO, and public company experience; M&A transaction support and international tax structures; SOX, financial infrastructure redesign and systems implementation. Result-driven and team-oriented professional with 15 years of combined experience in finance, accounting and financial reporting, and audit and assurance.',
                                linkedin: 'lbankurova'
                            })}

                            {rednderCard(props = {
                                name: 'Alex Paramonov',
                                id: 'alex-paramonov',
                                position: 'Software Architect',
                                avatar: '/docusaurus_img/team/ap.jpeg',
                                description: `I've wrote my first program in Basic when I was 6, and since then I'm completely fascinated with computers.
                                I really love all the electronic and digital devices, and my master's degree is engeneer of radio-electronic.
                                When I first took a programming side-job when I was 15, I've noticed how cool it is, and it's not even a job, since I was really happy doing this.
                                These days I work in Datagrok and I'm responsible for the architecture and Datagrok backend features, making applied developers and data-scientists experience smooth and productive.
                                My life-credo is to make people's life easier, and leave things better than they were before me.
                                When I'm not working, I usually make various electronic devices, compose music or play the electric guitar.

                                I really like to make things faster. When I worked on commercial banks, we used to make tools that completely automate processes usually taking days to seconds. Also, it's always exciting to participate in real-world problem solution projects, such as autism spectrum disorder research. I'm really proud that I took part in data processing pipeline development for this project.

                                I think, working in Datagrok requires you to really love computer science and have a holistic understanding of how computers work, together with the wish to make the world a better and more understandable place.`,
                                linkedin: 'alexaprm',
                                github: 'alex-aprm'
                            })}

                            {rednderCard(props = {
                                name: 'Sofiia Podolskaia',
                                id: 'sofiia-podolskaia',
                                position: 'DevOps Architect, Scrum Master',
                                avatar: '/docusaurus_img/team/sp.jpeg',
                                linkedin: 'sofiia-podolskaia',
                            })}

                            {rednderCard(props = {
                                name: 'Leonid Stolbov',
                                id: 'leonid-stolbov',
                                position: 'Research Scientist',
                                avatar: '/docusaurus_img/team/ls.jpeg',
                                linkedin: 'leonid-s-436742216',
                                github: 'StLeonidas'
                            })}

                            {rednderCard(props = {
                                name: 'Diana Onufriienko',
                                id: 'diana-onufriienko',
                                position: 'Software Engineer',
                                avatar: '/docusaurus_img/team/od.jpeg',
                                linkedin: 'diana-onufriienko',
                                github: 'onuf'
                            })}

                            {rednderCard(props = {
                                name: 'Konstantin Doncov',
                                id: 'konstantin-doncov',
                                position: 'Data Scientist',
                                avatar: '/docusaurus_img/team/kd.jpeg',
                                linkedin: 'konstantin-doncov-252064159',
                                github: 'konstantin-doncov'
                            })}

                            {rednderCard(props = {
                                name: 'Oleksii Sakhniuk',
                                id: 'oleksii-sakhniuk',
                                position: 'Quality Assurance Engineer',
                                avatar: '/docusaurus_img/team/os.jpeg',
                                linkedin: 'osakhniuk',
                                github: 'osakhniuk'
                            })}

                            {rednderCard(props = {
                                name: 'Aleksandr Tanas',
                                id: 'aleksandr-tanas',
                                position: 'Software Developer',
                                avatar: '/docusaurus_img/team/at.jpeg',
                                description:`Programming was my hobby back in my school days. My first degree was as an industrial electronics engineer. When graduated I worked as a professional software developer for about 8 years. Then I fall for genetics and studied bioinformatics at Yandex data analysis school simaltenoisly working on my PhD thesis in the field of molecular genetics DNA methylation screening methods. After the defence, I work at Research Centre for Medical Genetics studying mainly abnormal DNA methylation markers in cancer, maintaining the IT infrastructure of the laboratory, and developing scientific as well as management software. But war ruined everything and now I search for a job preferably as a software engineer.`,
                                linkedin: 'aleksandr-tanas',
                                github: 'tanas80'
                            })}

                            {rednderCard(props = {
                                name: 'Dmitry Illarionov',
                                id: 'dmitry-illarionov',
                                position: 'User Experience Engineer',
                                avatar: '/docusaurus_img/team/di.jpeg',
                                linkedin: 'illdm',
                                github: 'illarionow'
                            })}

                            {rednderCard(props = {
                                name: 'Kostiantyn Rudik',
                                id: 'kstiantyn-rudik',
                                position: 'Talent Manager',
                                avatar: '/docusaurus_img/team/kr.jpeg',
                                linkedin: '%D0%BA%D0%BE%D0%BD%D1%81%D1%82%D0%B0%D0%BD%D1%82%D0%B8%D0%BD-%D1%80%D1%83%D0%B4%D0%B8%D0%BA/',
                            })}

                            {rednderCard(props = {
                                name: 'Elena Zhuravskaya',
                                id: 'elena-zhuravskaya',
                                position: 'Junior Recruiter',
                                avatar: '/docusaurus_img/team/ez.jpeg',
                                linkedin: 'elena-zhuravskaya-2b2bb014a',
                            })}

                            {rednderCard(props = {
                                name: 'Volodymyr Dyma',
                                id: 'volodymyr-dyma',
                                position: 'Junior Data Scientist',
                                avatar: '/docusaurus_img/team/vd.jpeg',
                                linkedin: 'vdyma',
                                github: 'vdyma'
                            })}

                            {rednderCard(props = {
                                name: 'Denis Kryvchuk',
                                id: 'denis-kryvchuk',
                                position: 'Junior Developer',
                                avatar: '/docusaurus_img/team/dk.jpeg',
                                linkedin: 'denis-kryvchuk-2959801b5',
                                github: 'dionissqq'
                            })}

                            {rednderCard(props = {
                                name: 'Oleksandra Serhiienko',
                                id: 'oleksandra-serhiienko',
                                position: 'Junior Developer',
                                avatar: '/docusaurus_img/team/ss.jpeg',
                                linkedin: '%D0%B0%D0%BB%D0%B5%D0%BA%D1%81%D0%B0%D0%BD%D0%B4%D1%80%D0%B0-%D1%81%D0%B5%D1%80%D0%B3%D0%B8%D0%B5%D0%BD%D0%BA%D0%BE-674ab6239/',
                                github: 'Aleksashka11'
                            })}
                        </div>
                    </div>
                </section>
                {contactus()}
                {footer()}
                <a href="https://razomforukraine.org/" class="d4-link-external" style={{
                    'position': 'fixed',
                    'bottom': '0',
                    'left': '-80px',
                    'height': '84px',
                    'width': '300px',
                    'transform': 'rotate(45deg)',
                    'background': 'linear-gradient(-180deg,#1796FC 50%,#ffd500 0)',
                    'z-index': '99999',
                }}></a>
        </div >
    );
}
