#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Pmw',
    'forcefield',
    'logging',
    'molssi_util',
    'molssi_workflow',
    'pprint',
]
# 'math',
# 'random',

setup_requirements = [
    'pytest-runner',
    # TODO(paulsaxe): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'pytest',
    # TODO: put package test requirements here
]

setup(
    name='lammps_step',
    version='0.1.0',
    description="The LAMMPS step for a MolSSI workflow",
    long_description=readme + '\n\n' + history,
    author="Paul Saxe",
    author_email='psaxe@molssi.org',
    url='https://github.com/paulsaxe/lammps_step',
    packages=find_packages(include=['lammps_step']),
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='lammps_step',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
    entry_points={
        'org.molssi.workflow': [
            'LAMMPS = lammps_step:LAMMPSStep',
        ],
        'org.molssi.workflow.tk': [
            'LAMMPS = lammps_step:LAMMPSStep',
        ],
        'org.molssi.workflow.lammps': [
            'Custom = lammps_step:CustomStep',
            'Energy = lammps_step:EnergyStep',
            'Initialization = lammps_step:InitializationStep',
            'Minimization = lammps_step:MinimizationStep',
            'NVE = lammps_step:NVEStep',
            'NVT = lammps_step:NVTStep',
            'Velocities = lammps_step:VelocitiesStep',
        ],
        'org.molssi.workflow.lammps.tk': [
            'Custom = lammps_step:CustomStep',
            'Energy = lammps_step:EnergyStep',
            'Initialization = lammps_step:InitializationStep',
            'Minimization = lammps_step:MinimizationStep',
            'NVE = lammps_step:NVEStep',
            'NVT = lammps_step:NVTStep',
            'Velocities = lammps_step:VelocitiesStep',
        ],
    }
)
