from setuptools import setup


def readme():
	with open('README.rst') as f:
		return f.read()


setup(name='tempertonlab_utils',
      version='0.1',
      description='Functions frequently used by the Temperton Lab',
      long_description=readme(),
      classifiers=[
	      'Development Status :: 3 - Alpha',
	      'License :: OSI Approved :: MIT License',
	      'Programming Language :: Python :: 3.7',
	      'Topic :: Genomics',
      ],
      url='https://github.com/btemperton/tempertonlab_utils',
      author='Ben Temperton',
      author_email='btemperton@gmail.com',
      license='MIT',
      packages=['tempertonlab_utils'],
      install_requires=[
	      'numpy',
	      'pandas',
	      'kpal',
	      'scikit-learn',
	      'biopython'],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['bin/rename_fastas']
      )
