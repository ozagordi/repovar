#!/usr/bin/env python

import setuptools

setuptools.setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    install_requires=['setuptools_scm'],
    name='repovar',
      description='Report minority variants found in genetically heterogeneous samples',
      url='https://github.com/ozagordi/repovar',
      author='Osvaldo Zagordi',
      author_email='firstname.lastname@gmail.com',
      packages = setuptools.find_packages(),
      scripts = ['bin/repovar'],
      data_files = [('repovar/templates',
        ['repovar/templates/template.tex'])],
      entry_points={
          'console_scripts': [
              'repovar = repovar.__main__:main'
          ]
    }
)
