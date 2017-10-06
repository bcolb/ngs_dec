import glob
import os
import subprocess

from distutils.core import Command
from setuptools import setup, find_packages


class CheckVersion(Command):
    description = 'Confirm that the stored package version is correct'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        with open('ngs_dec/data/ver') as f:
            stored_version = f.read().strip()

        git_version = subprocess.check_output(
            ['git', 'describe', '--tags', '--dirty']).strip()

        assert stored_version == git_version
        print 'the current version is', stored_version

subprocess.call(
    ('mkdir -p mkvenv/data && '
     'git describe --tags --dirty > ngs_dec/data/ver.tmp '
     '&& mv ngs_dec/data/ver.tmp ngs_dec/data/ver '
     '|| rm -f ngs_dec/data/ver.tmp'),
    shell=True, stderr=open(os.devnull, "w"))

from ngs_dec import __version__
package_data = glob.glob('data/*')

params = {'author': 'Your name',
          'author_email': 'Your email',
          'description': 'Package description',
          'name': 'ngs_dec',
          'packages': find_packages(),
          'package_dir': {'ngs_dec': 'ngs_dec'},
          'entry_points': {
              'console_scripts': ['ngs_dec = ngs_dec.scripts.main:main']
          },
          'version': __version__,
          'setup_requires': 'numpy==1.13.1',
          'install_requires': [
              'pandas==0.20.3',
              'numpy==1.13.1',
              'scipy==0.19.1'],
          'package_data': {'ngs_dec': package_data},
          'test_suite': 'tests',
          'cmdclass': {'check_version': CheckVersion}
          }

setup(**params)
