# Systematic Review PDF Extractor

A web application that helps researchers extract information from PDF files for systematic reviews.

## Features

- Upload multiple PDF files
- Upload a protocol file
- Extract key information from PDFs
- Generate an Excel file with extracted data
- Download results directly from the browser

## How to Use

1. Upload your PDF files (multiple files supported)
2. Upload your protocol file
3. Click "Start Extraction"
4. Wait for the processing to complete
5. Download the Excel file with the results

## Local Development

To run the app locally:

1. Install the requirements:
```bash
pip install -r requirements.txt
```

2. Run the app:
```bash
streamlit run src/agent/pdf_extractor_streamlit.py
```

## Online Version

The app is available online at: [Your Streamlit Cloud URL will be here]

## Requirements

- Python 3.8+
- Streamlit
- Pandas
- PyPDF2
- Openpyxl
