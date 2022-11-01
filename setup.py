import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snapper-ont",
    version="0.3.1",
    author="D.N. Konanov",
    author_email="konanovdmitriy@gmail.com",
    description="Nanopore-based methylation sites caller",
    long_description="snapper",
    long_description_content_type="",
    url="https://github.com/DNKonanov/Snapper",
    project_urls={
        "Bug Tracker": "https://github.com/DNKonanov/Snapper",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    packages=['snapper', 'snapper.src'],
    install_requires=[
        'h5py',
        'biopython',
        'matplotlib',
        'scipy',
        'seaborn'
    ],
    entry_points={
        'console_scripts': [
            'snapper=snapper.snapper:main'
        ]
    }
)
