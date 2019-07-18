#!/usr/bin/env python

from distutils.core import setup
import os

def fread_db_library(directory):

    paths = []
    for path, _, templates in os.walk(directory):
        for template in templates:
            template = os.path.join(path, template).split("lib/python/saab_plus/")[-1]
            paths.append(template)
    return paths

setup(name='saab_plus',
      version='1.0a',
      description='Structural Diversity of B-Cell Repertoire along the B-cell Differentiation axis in Humans',
      author='Aleksandr Kovaltsuk',
      author_email='opig@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk//resources',

      packages=["saab_plus",
                "saab_plus.aboss_utils",
                "saab_plus.code",
                "saab_plus.code.Common",
                "saab_plus.code.DataManagement",
                "saab_plus.code.Alignment",
                "saab_plus.code.Alignment.FREAD",
                "saab_plus.code.Alignment.FREAD.tools",
                "saab_plus.code.Alignment.FREAD.prosci",
                "saab_plus.code.Alignment.FREAD.prosci.loops",
                "saab_plus.code.Alignment.FREAD.prosci.util"],
      package_dir={"saab_plus": "lib/python/saab_plus",
                   "saab_plus.aboss_utils": "lib/python/saab_plus/aboss_utils",
                   "saab_plus.code": "lib/python/saab_plus/code",
                   "saab_plus.code.Common": "lib/python/saab_plus/code/Common",
                   "saab_plus.code.DataManagement": "lib/python/saab_plus/code/DataManagement",
                   "saab_plus.code.Alignment": "lib/python/saab_plus/code/Alignment",
                   "saab_plus.code.Alignment.FREAD": "lib/python/saab_plus/code/Alignment/FREAD",
                   "saab_plus.code.Alignment.FREAD.tools": "lib/python/saab_plus/code/Alignment/FREAD/tools",
                   "saab_plus.code.Alignment.FREAD.prosci": "lib/python/saab_plus/code/Alignment/FREAD/prosci",
                   "saab_plus.code.Alignment.FREAD.prosci.loops": "lib/python/saab_plus/code/Alignment/FREAD/prosci/loops",
                   "saab_plus.code.Alignment.FREAD.prosci.util": "lib/python/saab_plus/code/Alignment/FREAD/prosci/util"},
      package_data = {"saab_plus": fread_db_library("lib/python/saab_plus/data")+ 
                                   ["data_esst/esst.txt"]
                                    },
      scripts=["bin/SAAB_PLUS", 
               "bin/SAAB_PLUS_DIAG"],
      license="MIT",
     )


import sys
if sys.argv[1] != "install":
    sys.exit(0)


