import sys
import django
from os import path as osp


def rel_path(*p): return osp.normpath(osp.join(rel_path.path, *p))


rel_path.path = osp.abspath(osp.dirname(__file__))
this = osp.splitext(osp.basename(__file__))[0]

from django.conf import settings

SETTINGS = dict(
    SITE_ID=1,
    DATABASES={},
    DEBUG=True,
    TEMPLATE_DEBUG=True,
    ROOT_URLCONF=this
)
SETTINGS['TEMPLATE_DIRS'] = (rel_path(),),
SETTINGS['DATABASES'] = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': rel_path('db')
    }
}

SETTINGS['INSTALLED_APPS'] = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    # 'django.contrib.staticfiles',
    'django.contrib.admin',
    'corsheaders',
    'rest_framework'
)

if not settings.configured:
    settings.configure(**SETTINGS)

django.setup()

from django.db import models
from django.contrib import admin

class Smiles(models.Model):
    smiles = models.TextField()
    numerical_data = models.TextField()
    class Meta:
        app_label = this
    __module__ = this

try:
    admin.site.register(Smiles)
except admin.sites.AlreadyRegistered:
    pass

admin.autodiscover()

from rest_framework import serializers

class SmilesSerializer (serializers.ModelSerializer):
    class Meta:
        model = Smiles
        fields = ('smiles', 'numerical_data')

from django.core.urlresolvers import reverse
from rest_framework import viewsets, status
from rest_framework.decorators import list_route
from rest_framework.response import Response
from rest_framework.settings import api_settings
from rest_framework_csv.parsers import CSVParser
from rest_framework_csv.renderers import CSVRenderer

class SmilesViewSet(viewsets.ModelViewSet):
    queryset = Smiles.objects.all()
    parser_classes = (CSVParser,) + tuple(api_settings.DEFAULT_PARSER_CLASSES)
    renderer_classes = (CSVRenderer,) + tuple(api_settings.DEFAULT_RENDERER_CLASSES)
    serializer_class = SmilesSerializer

    def get_renderer_context(self):
        context = super(SmilesViewSet, self).get_renderer_context()
        context['header'] = (
            self.request.GET['fields'].split(',')
            if 'fields' in self.request.GET else None)
        return context

    @list_route(methods=['POST'])
    def df_upload(self, request, *args, **kwargs):
        return Response(handle_uploaded_file(request.data, request.query_params.get('models')))

import sklearn.externals.joblib
import numpy as np
import pandas as pd
import csv
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from pychem import pychem
from pychem.pychem import Chem
from pychem import constitution
from pychem import topology
from pychem import connectivity as con
from pychem import kappa
from pychem import estate
from pychem import basak
from pychem import moran, geary
from pychem import molproperty as mp
from pychem import charge
from pychem import moe

"""Dictionary that defines number of bits for each model or list of needed descriptors"""
dict_bits_desc = {
    "Pgp-Inhibitor/Pgp-Inhibitor": 2048,
    "Pgp-Substrate/Pgp-Substrate": 2048,
    "HIA": 167,
    "F(20%)": 167,
    "F(30%)": 2048,
    "Ames": 167,
    "SkinSen": 167,
    "BBB/BBB": 2048,
    "CYP1A2-Inhibitor/CYP1A2-Inhibitor": 2048,
    "CYP1A2-Substrate": 1024,
    "CYP3A4-Inhibitor/CYP3A4-Inhibitor": 2048,
    "CYP3A4-Substrate": 1024,
    "CYP2C19-Inhibitor/CYP2C19-Inhibitor": 2048,
    "CYP2C19-Substrate": 1024,
    "CYP2C9-Inhibitor/CYP2C9-Inhibitor": 2048,
    "CYP2C9-Substrate": 1024,
    "CYP2D6-Inhibitor": 1024,
    "CYP2D6-Substrate": 1024,
    "Clearance/Clearance": ['nsulph', 'VSAEstate8', 'QNmin', 'IDET', 'ndb', 'slogPVSA2', 'MATSv5', 'S32', 'QCss',
                 'bcutm4', 'S9', 'bcutp8', 'Tnc', 'nsb', 'Geto', 'bcutp11', 'S7', 'MATSm2', 'GMTIV',
                 'nhet', 'MATSe5', 'CIC0', 'bcutp3', 'Gravto', 'EstateVSA9', 'MATSe3', 'MATSe5', 'UI',
                 'S53', 'J', 'bcute1', 'MRVSA9', 'PEOEVSA0', 'MATSv2', 'IDE', 'AWeight', 'IC0', 'S16',
                 'bcutp1', 'PEOEVSA12'],
    "T/T": ['MATSv5', 'Gravto', 'Chiv3c', 'PEOEVSA7', 'knotp', 'bcutp3', 'bcutm9', 'EstateVSA3',
                'MATSp1', 'bcutp11', 'VSAEstate7', 'IC0', 'UI', 'Geto', 'QOmin', 'CIC0', 'dchi3',
                'MATSp4', 'bcutm4', 'Hatov', 'MATSe4', 'CIC6', 'Chiv4', 'EstateVSA9', 'MATSv2', 'nring',
                'bcute1', 'VSAEstate8', 'MRVSA9', 'PEOEVSA6', 'SIC1', 'bcutp8', 'MATSp6', 'QCss', 'J',
                'IDE', 'CIC2', 'Hy', 'MRVSA6', 'naro', 'SPP', 'EstateVSA7', 'bcutv10', 'S12', 'LogP2',
                'bcutp2', 'CIC3', 'S17', 'LogP', 'bcutp1'],
    "hERG": ['ndb', 'nsb', 'ncarb', 'nsulph', 'naro', 'ndonr', 'nhev', 'naccr', 'nta', 'nring',
                   'PC6', 'GMTIV', 'AW', 'Geto', 'BertzCT', 'J', 'MZM2', 'phi', 'kappa2', 'MATSv1',
                   'MATSv5', 'MATSe4', 'MATSe5', 'MATSe6', 'TPSA', 'Hy', 'LogP', 'LogP2', 'UI', 'QOss',
                   'SPP', 'LDI', 'Qass', 'QOmin', 'QNmax', 'Qmin', 'Mnc', 'EstateVSA7', 'EstateVSA0',
                   'EstateVSA3', 'PEOEVSA0', 'PEOEVSA6', 'MRVSA5', 'MRVSA4', 'MRVSA3', 'MRVSA6',
                   'slogPVSA1'],
    "H-HT": ['ndb', 'nsb', 'nnitro', 'naro', 'ndonr', 'nhet', 'nhev', 'nring', 'PC6', 'GMTIV',
                  'Geto', 'IDE', 'Arto', 'Hatov', 'BertzCT', 'Getov', 'J', 'MZM2', 'MZM1', 'phi',
                  'kappa3', 'kappa2', 'kappam3', 'MATSp4', 'MATSp6', 'MATSv3', 'MATSv2', 'MATSv5',
                  'MATSv7', 'MATSv6', 'MATSm4', 'MATSm5', 'MATSm6', 'MATSm2', 'MATSm3', 'MATSe4',
                  'MATSe5', 'MATSe6', 'MATSe1', 'MATSe2', 'MATSe3', 'MATSp3', 'MATSp2', 'TPSA', 'LogP',
                  'LogP2', 'UI', 'QNmin', 'QOss', 'QHss', 'SPP', 'LDI', 'Qass', 'QCmax', 'QOmax', 'Tpc',
                  'QOmin', 'QCss', 'QHmax', 'Rnc', 'Rpc', 'Qmin', 'Mnc', 'EstateVSA9', 'EstateVSA4',
                  'EstateVSA5', 'EstateVSA6', 'EstateVSA7', 'EstateVSA0', 'EstateVSA1', 'EstateVSA2',
                  'EstateVSA3', 'PEOEVSA11', 'PEOEVSA2', 'PEOEVSA1', 'PEOEVSA7', 'PEOEVSA6', 'PEOEVSA5',
                  'MRVSA5', 'MRVSA4', 'PEOEVSA9', 'PEOEVSA8', 'MRVSA3', 'MRVSA2', 'MRVSA9', 'MRVSA6',
                  'slogPVSA2', 'slogPVSA4', 'slogPVSA5'],
    "LD50": ['ATSm1', 'ATSm2', 'ATSm3', 'ATSm4', 'ATSm6', 'AWeight', 'Chi4c', 'Chiv3', 'Chiv4',
                 'Chiv4c', 'Chiv4pc', 'DS', 'Gravto', 'IC0', 'IC1', 'MRVSA9', 'QCmax', 'QNss', 'QOmin',
                 'Qmax', 'S46', 'Smax45', 'Smin', 'Smin45', 'VSAEstate7', 'Weight', 'bcutm1', 'bcutm2',
                 'bcutp1', 'nhet', 'nphos', 'slogPVSA11'],
    "PPB": ['ncarb', 'naro', 'GMTIV', 'AW', 'Geto', 'Arto', 'BertzCT', 'J', 'kappam3', 'MATSv1',
                  'MATSe1', 'Hy', 'LogP', 'LogP2', 'UI', 'QHss', 'LDI', 'QNss', 'Rpc', 'Mnc', 'PEOEVSA10',
                  'PEOEVSA0', 'PEOEVSA6', 'PEOEVSA5', 'PEOEVSA4', 'slogPVSA10', 'MRVSA6', 'slogPVSA0',
                  'slogPVSA1', 'slogPVSA5'],
    "VD/VD": ['GMTIV', 'UI', 'MATSe1', 'MATSp1', 'Chiv4', 'MATSm2', 'S12', 'dchi3', 'IDE', 'PEOEVSA7',
                 'bcutp1', 'bcutm9', 'SIC1', 'MRVSA6', 'IC1', 'QNmax', 'CIC0', 'PEOEVSA6', 'MATSe4',
                 'VSAEstate8', 'Geto', 'EstateVSA3', 'MRVSA5', 'LogP2', 'Tnc', 'S7', 'SPP', 'QOmin',
                 'EstateVSA7', 'LogP', 'QNmin', 'MRVSA9', 'S19', 'MATSv2', 'nsulph', 'S17', 'S9', 'ndb',
                 'AWeight', 'QCss', 'EstateVSA9', 'Hy', 'S16', 'IC0', 'S30']
}

def calculate_descriptors(smile):
    """Calculate descriptors for the input smile

    :param smile: input smile
    :return: dict with calculated descriptors (key=descriptor name, value=calculated descriptor value)
    """
    res_dict = {}
    mol = Chem.MolFromSmiles(smile)
    res_dict['MolWeight'] = constitution.CalculateMolWeight(mol)
    res_dict['HeavyAtomNumber'] = constitution.CalculateHeavyAtomNumber(mol)
    res_dict['Path2'] = constitution.CalculatePath2(mol)
    res_dict['Constitutional'] = constitution.GetConstitutional(mol)
    res_dict['Balaban'] = topology.CalculateBalaban(mol)
    res_dict['MZagreb1'] = topology.CalculateMZagreb1(mol)
    res_dict['Harary'] = topology.CalculateHarary(mol)
    res_dict['Topology'] = topology.GetTopology(mol)
    res_dict['Chi2'] = con.CalculateChi2(mol)
    res_dict['MeanRandic'] = con.CalculateMeanRandic(mol)
    res_dict['Connectivity'] = con.GetConnectivity(mol)
    res_dict['Kappa1'] = kappa.CalculateKappa1(mol)
    res_dict['Kappa2'] = kappa.CalculateKappa2(mol)
    res_dict['GetKappa'] = kappa.GetKappa(mol)
    drug = pychem.PyChem2d()
    drug.ReadMolFromSmile(smile)
    res_dict['Bcut'] = drug.GetBcut()
    res_dict['HeavyAtomEstate'] = estate.CalculateHeavyAtomEState(mol)
    res_dict['MaxEstate'] = estate.CalculateMaxEState(mol)
    res_dict['HalogenEstate'] = estate.CalculateHalogenEState(mol)
    res_dict['GetEstate'] = estate.GetEstate(mol)
    res_dict['BasakCIC1'] = basak.CalculateBasakCIC1(mol)
    res_dict['BasakSIC1'] = basak.CalculateBasakSIC1(mol)
    res_dict['BasakSIC3'] = basak.CalculateBasakSIC3(mol)
    res_dict['GetBasak'] = basak.Getbasak(mol)
    res_dict['MoranAutoVolume'] = moran.CalculateMoranAutoVolume(mol)
    res_dict['GeatyAutoMass'] = geary.CalculateGearyAutoMass(mol)
    res_dict['MoranAuto'] = moran.GetMoranAuto(mol)
    res_dict['logP'] = mp.CalculateMolLogP(mol)
    res_dict['MR'] = mp.CalculateMolMR(mol)
    res_dict['TPSA'] = mp.CalculateTPSA(mol)
    res_dict['MolProperty'] = mp.GetMolecularProperty(mol)
    res_dict['DipoleIndex'] = charge.CalculateLocalDipoleIndex(mol)
    res_dict['SquareCharge'] = charge.CalculateAllSumSquareCharge(mol)
    res_dict['GetCharge'] = charge.GetCharge(mol)
    res_dict['moeTPSA'] = moe.CalculateTPSA(mol)
    res_dict['PEOEVSA'] = moe.CalculatePEOEVSA(mol)
    res_dict['GetMoe'] = moe.GetMOE(mol)
    for k, v in res_dict.items():
        if type(v) == dict:
            for k1, v1 in v.items():
                res_dict[k1] = v1
            del res_dict[k]
    return res_dict

def get_descriptor_vector(smiles, desc_names):
    """Calculate descriptor vectors for some excretion, toxicity and distribution models

    :param smiles: list of smiles
    :param desc_names: list of descriptor names that need to be calculated
    :return: list with calculated descriptors for each smile
    """
    desc_vector = []
    for smile in smiles:
        desc_vector_smile = []
        desc_dict = calculate_descriptors(smile)
        for name in desc_names:
            if name in desc_dict.keys():
                desc_vector_smile.append(desc_dict[name])
        desc_vector.append(desc_vector_smile)
    return desc_vector

def getMACCS(smiles):
    """Calculate MACCS fingerprints for the list of smiles

    :param smiles: list of smiles
    :return: list with calculated MACCS fingerprints (as a bit vector) for each smile
    """
    result = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        fpM = MACCSkeys.GenMACCSKeys(mol)
        bs2 = ','.join(fpM.ToBitString())
        result.append(bs2.split(','))
    return result

def getECFP(smiles, radius, nBits):
    """Calculate ECFP fingerprints depending on the radius and size of the bitset

    :param smiles: list of smiles
    :param radius: radius
    :param nBits: number of bits to be returned
    :return: list with calculated ECFP fingerprints (as a bit vector) for each smile
    """
    result = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        fpECFP = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useChirality=False)
        result.append(fpECFP)
    return result

def handle_uploaded_file(f, models):
    """Calculate ADMET properties

    :param f: input file content
    :param models: list of models to be used for calculation
    :return: zip(predict_label, predict_proba)
    """
    smiles = []
    encoded_smiles = []
    result = []
    for row in f:
        smiles.append(list(row.values()))
    for smile in smiles:
        encoded_smiles.append(smile[0].encode('utf8'))
    current_path = os.path.split(os.path.realpath(__file__))[0]
    for model in models.split(","):
        cf = sklearn.externals.joblib.load(current_path + '/' + model)
        fingerprint_content = lambda bits: getMACCS(encoded_smiles) if bits == 167 \
                              else (getECFP(encoded_smiles, 1, 2048) if bits == 2048 \
                              else (getECFP(encoded_smiles, 2, 1024) if bits == 1024 \
                              else get_descriptor_vector(encoded_smiles, bits)))
        des_list = np.array(fingerprint_content(dict_bits_desc[str(model)[:-4]]))
        y_predict_label = cf.predict(des_list)
        y_predict_proba = cf.predict_proba(des_list)
        result.append(zip(y_predict_label, y_predict_proba))
    return result

from django.conf.urls import url, include
from django.contrib import admin
from rest_framework import routers

router = routers.DefaultRouter()
router.register(r'smiles', SmilesViewSet)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^admin/', admin.site.urls),
]

# CLI

if __name__ == "__main__":
    # make this script runnable like a normal `manage.py` command line script.
    from django.core import management
    management.execute_from_command_line()