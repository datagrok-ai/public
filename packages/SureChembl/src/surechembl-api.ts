// SureChemBL API Types

export interface SureChemblBatchResponse {
    status: string;
    timestamp: string;
    error_message: string;
    data: PatentDocument[];
}

export interface PatentDocument {
    doc_id: string;
    doc_version: string;
    date_added: string;
    contents: {
        patentDocument: PatentContent;
    };
}

export interface PatentContent {
    bibliographicData: BibliographicData;
    family: string[];
    amendedClaims: string;
    amendedClaimsStatement: string;
    searchReportData: string;
    legalStatus: LegalStatus[];
    revisionHistory: string;
    copyright: {
        t_sys_data: string;
    };
    pdfUrl: string;
    abstracts: Abstract[];
    descriptions: Description[];
    claimResponses: ClaimResponse[];
    published: string;
    drawingResponse: DrawingResponse;
    hasPDF: boolean;
}

export interface BibliographicData {
    publicationReference: PublicationReference[];
    usSlrFlag: string;
    applicationReference: ApplicationReference[];
    priorityClaims: PriorityClaim[];
    datesOfPublicAvailability: string;
    rule47Flag: string;
    termOfGrant: TermOfGrant[];
    technicalData: TechnicalData;
    relatedDocuments: RelatedDocument[];
    parties: Party[];
    internationalConventionData: string;
    designatedStates: string;
    officeSpecificData: string;
}

export interface PublicationReference {
    fvid: string;
    ucid: string;
    entityStatus: string;
    documentId: DocumentId[];
}

export interface DocumentId {
    country: {
        content: string;
    };
    docNumber: string;
    kind: string;
    date: string;
    lang: string;
}

export interface ApplicationReference {
    ucid: string;
    isRepresentative: string;
    usArtUnit: string;
    usSeriesCode: string;
    documentIdList: DocumentId[];
}

export interface PriorityClaim {
    priorityClaims: PriorityClaimDetail[];
}

export interface PriorityClaimDetail {
    mxwId: string;
    ucid: string;
    loadSource: string;
    documentId: DocumentId;
    kind: string;
}

export interface TermOfGrant {
    'us-term-extension': {
        $: string;
    };
}

export interface TechnicalData {
    classificationsIpc: string;
    classificationLocarno: string;
    classificationsIpcr: ClassificationIpcr[];
    classificationNational: ClassificationNational[];
    classificationEcla: ClassificationEcla[];
    classificationsCpc: ClassificationCpc[];
    fTermInfo: string;
    citations: Citation[];
    figures: Figure[];
    fieldOfSearches: FieldOfSearch[];
    inventionTitles: InventionTitle[];
}

export interface ClassificationIpcr {
    classificationIpcr: {
        classification: string;
        mxwId: string;
        loadSource: string;
    }[];
}

export interface ClassificationNational {
    country: string;
    mainClassification: {
        content: string;
        mxwId: string;
        loadSource: string;
    };
    'further-classification': FurtherClassification[];
    furtherClassifications: FurtherClassification[];
}

export interface FurtherClassification {
    mxwId: string;
    loadSource: string;
    value: string;
}

export interface ClassificationEcla {
    classificationSymbol: {
        mxwId: string;
        loadSource: string;
        scheme: string;
        classification: string;
    }[];
}

export interface ClassificationCpc {
    classificationCpc: {
        scheme: string;
        classification: string;
        mxwId: string;
        loadSource: string;
    }[];
}

export interface Citation {
    patentCitations: {
        patcit: PatentCitation[];
    };
}

export interface PatentCitation {
    mxwId: string;
    loadSource: string;
    ucid: string;
    documentId: DocumentId;
    sources: {
        source: {
            name: string;
            createdByNpl: string;
        }[];
    };
}

export interface Figure {
    'number-of-drawing-sheets': {
        $: string;
    };
    'number-of-figures': {
        $: string;
    };
}

export interface FieldOfSearch {
    'classification-national': {
        country: string;
        'further-classification': string;
        text: string;
        mxwId: string;
        loadSource: string;
    }[];
}

export interface InventionTitle {
    lang: string;
    title: string;
    mxwid: string;
    annotation: {
        rank: number;
        field: string;
        chemicalAnnotations: ChemicalAnnotation[];
        diseaseAnnotations: any[];
        modeOfActionAnnotations: any[];
        targetAnnotations: any[];
    };
}

export interface ChemicalAnnotation {
    name: string;
    start: number;
    end: number;
    category: string;
    globalFrequency: number;
    chemicalIds: number[];
}

export interface RelatedDocument {
    relation: {
        '@type': string;
        'document-id': DocumentId;
    };
}

export interface Party {
    applicants?: {
        applicant: Applicant[];
    };
    inventors?: {
        inventor: Inventor[];
    };
    assignees?: {
        assignee: Assignee[];
    };
    agents?: {
        agent: Agent[];
    };
    assigneeHistory?: {
        reassignments: {
            reassignment: Reassignment[];
        };
    };
    examiners?: {
        examiner: Examiner[];
    };
}

export interface Applicant {
    mxwId: string;
    loadSource: string;
    sequence: string;
    format: string;
    addressbook: {
        lastName?: string;
        address?: {
            country?: string;
        };
    };
}

export interface Inventor {
    mxwId: string;
    loadSource: string;
    sequence: string;
    format: string;
    addressbook: {
        lastName?: string;
        firstName?: string;
        address?: {
            city?: string;
            state?: string;
            country?: string;
        };
    };
}

export interface Assignee {
    mxwId: string;
    loadSource: string;
    sequence: string;
    format: string;
    usAssigneeType?: string;
    addressbook: {
        name?: string;
        lastName?: string;
        address?: {
            address1?: string;
            address2?: string;
            address3?: string;
            city?: string;
            state?: string;
            postcode?: string;
            country?: string;
        };
    };
}

export interface Agent {
    mxwId: string;
    loadSource: string;
    sequence: string;
    format: string;
    addressbook: {
        name?: string;
    };
}

export interface Reassignment {
    mxwId: string;
    loadSource: string;
    reelFrame: string;
    recordedDate: string;
    assignees: {
        assignee: Assignee[];
    };
    assignors: {
        assignor: Assignor[];
    };
    correspondenceAddress: {
        mxwId: string;
        loadSource: string;
        format: string;
        addressbook: {
            lastName: string;
            address: {
                address1?: string;
                address2?: string;
                address3?: string;
            };
        };
    };
    conveyance: string;
}

export interface Assignor {
    mxwId: string;
    loadSource: string;
    format: string;
    executionDate: string;
    addressbook: {
        lastName: string;
    };
}

export interface Examiner {
    examinerLevel: string;
    mxwId: string;
    loadSource: string;
    sequence: string;
    format: string;
    addressbook: {
        lastName: string;
        department: string;
        firstName: string;
    };
}

export interface LegalStatus {
    'legal-event': LegalEvent[];
}

export interface LegalEvent {
    '@country': string;
    '@code': string;
    '@date': string;
    '@mxw-id': string;
    '@impact'?: string;
    'legal-event-body': {
        'event-title': {
            $: string;
        };
        'event-attributes': {
            'event-attribute': {
                'event-attribute-label': {
                    $: string;
                }[];
                'event-attribute-value': {
                    $: string;
                }[];
            };
        };
    };
}

export interface Abstract {
    mxwid: string;
    lang: string;
    section: PatentSection;
}

export interface Description {
    mxwid: string;
    lang: string;
    loadSource: string;
    section: PatentSection;
}

export interface ClaimResponse {
    lang: string;
    loadSource: string;
    mxwid: string;
    section: PatentSection;
}

export interface PatentSection {
    content: string;
        annotationProblem?: {
            annotationProblemsEnum: string;
            annotation: ChemicalAnnotation;
        };
    annotations: any[];
}

export interface DrawingResponse {
    figures: FigureDetail[];
    lang: string;
    id: string;
    mxwId: string;
    copyrightRef: string;
}

export interface FigureDetail {
    id: string;
    num: string;
    figureLabels: string;
    imageID: string;
    file: string;
    width: string;
    height: string;
    imgContent: string;
    imgFormat: string;
    orientation: string;
    inline: string;
    alt: string;
    url: string;
    path: string;
}
