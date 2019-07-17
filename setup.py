#!/usr/bin/env python

from distutils.core import setup
import distutils
import os

def fread_library(directory):
    paths = []
    for path, _, templates in os.walk(directory):
        for template in templates:
            paths.append(os.path.join(path, template))
    return paths

data_paths = fread_library("data/")

setup(name='saab_plus',
      version='1.0',
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
      package_dir={"saab_plus": ".",
                   "saab_plus.aboss_utils": "aboss_utils",
                   "saab_plus.code": "code",
                   "saab_plus.code.Common": "code/Common",
                   "saab_plus.code.DataManagement": "code/DataManagement",
                   "saab_plus.code.Alignment": "code/Alignment",
                   "saab_plus.code.Alignment.FREAD": "code/Alignment/FREAD",
                   "saab_plus.code.Alignment.FREAD.tools": "code/Alignment/FREAD/tools",
                   "saab_plus.code.Alignment.FREAD.prosci": "code/Alignment/FREAD/prosci",
                   "saab_plus.code.Alignment.FREAD.prosci.loops": "code/Alignment/FREAD/prosci/loops",
                   "saab_plus.code.Alignment.FREAD.prosci.util": "code/Alignment/FREAD/prosci/util"},
      package_data={'saab_plus': data_paths +["code/Alignment/FREAD/esst.txt"]
                    },
      license="MIT",
     )


import sys
if sys.argv[1] != "install":
    sys.exit(0)


