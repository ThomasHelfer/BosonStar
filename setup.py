from setuptools import setup
from setuptools import find_packages

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='bosonstar',
	version='0.1',
	description='HBA with GPR and PCA',
	long_description=readme(),
	classifiers=[
		'Development Status :: 3 - Alpha',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2.7',
		'Topic :: Gravitational Wave :: Physics',
	],
	keywords='Boson Stars',
#	url='http://github.com/storborg/funniest',
	author='ThomasHelfer',
	author_email='thomashelfer@live.de',
	license='MIT',
	packages=find_packages(exclude = ["tests"]),
	install_requires=[
		'numpy',
		'scipy',
                'matplotlib',
	],
#        python_requires='>=3.5 ',
	include_package_data=True,
	zip_safe=False)
