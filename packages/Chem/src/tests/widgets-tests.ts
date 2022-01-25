import {category, test} from "@datagrok-libraries/utils/src/test";
import {drugLikenessWidget} from '../widgets/drug-likeness';
import {identifiersWidget} from '../widgets/identifiers';
import {molfileWidget} from '../widgets/molfile';
import {propertiesWidget} from '../widgets/properties';
import {structuralAlertsWidget} from '../widgets/structural-alerts';
import {structure2dWidget} from '../widgets/structure2d';
import {structure3dWidget} from '../widgets/structure3d';
import {toxicityWidget} from '../widgets/toxicity';

category('Chem: Widgets', () => {
    const molStr = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';

    test('drug-likeness', async () => {
        drugLikenessWidget(molStr);
    });

    test('identifiers', async () => {
        identifiersWidget(molStr);
    });

    test('molfile', async () => {
        molfileWidget(molStr);
    });

    test('properties', async () => {
        propertiesWidget(molStr);
    });

    test('structural-alerts', async () => {
        structuralAlertsWidget(molStr);
    });

    test('structure-2d', async () => {
        structure2dWidget(molStr);
    });

    test('structure-3d', async () => {
        structure3dWidget(molStr);
    });

    test('toxicity', async () => {
        toxicityWidget(molStr);
    });
});
