import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requirements = [
        req.strip()
        for req in f.readlines()
        if not req.startswith("#") and req.__contains__("==")
    ]

setuptools.setup(
    name='core_poly_match',
    version='0.0.1',
    author='Huimei Feng',
    author_email='fenghm@nao.cas.cn',
    description='FGS guider matching tool of CSST',  # short description
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://csst-tb.bao.ac.cn/code/csst-utils/csst_fgs_match',
    project_urls={
        'Source': 'https://csst-tb.bao.ac.cn/code/csst-utils/csst_fgs_match',
    },
    packages=setuptools.find_packages(),
    license='MIT',
    classifiers=["Development Status :: 5 - Production/Stable",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: MIT License",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python :: 3.9",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Astronomy"],
    # package_dir={'core_poly_match': 'core_poly_match'},
    # include_package_data=True,
    include_package_data=False,
    python_requires='>=3.9',
    install_requires=requirements,
)


