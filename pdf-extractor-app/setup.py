from setuptools import setup, find_packages

setup(
    name="pdf-extractor-app",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "streamlit==1.45.1",
        "pandas==2.3.0",
        "PyPDF2==3.0.1",
        "openpyxl==3.1.5",
        "requests==2.31.0",
        "python-dotenv==1.0.0",
    ],
) 