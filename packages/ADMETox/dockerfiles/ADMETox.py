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

"""Dictionary that defines number of bits for each model"""
dict_bits = {
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
    "Clearance/Clearance": 40,
    "T/T": 50,
    "hERG": 47,
    "H-HT": 89,
    "LD50": 32,
    "PPB": 30,
    "VD/VD": 45
}

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
        cf = sklearn.externals.joblib.load(current_path + '/static/models/' + model)
        fingerprint_content = lambda bits: getMACCS(encoded_smiles) if bits == 167 \
                              else (getECFP(encoded_smiles, 1, 2048) if bits == 2048 \
                              else getECFP(encoded_smiles, 2, 1024))
        des_list = np.array(fingerprint_content(dict_bits[str(model)[:-4]]))
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