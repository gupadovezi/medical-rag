# AI-Powered PDF Extractor

This application helps you extract and analyze information from PDF files using AI. It uses OpenRouter's API to process research papers and extract structured information.

## Features

- Upload multiple PDF files
- Extract key information using AI
- Generate Excel file with structured data
- Analyze findings across multiple papers
- Identify common themes and research gaps

## Online Version

The app is deployed at: [Your Streamlit Cloud URL will appear here]

## Local Development

1. Clone the repository:
```bash
git clone https://github.com/gupadovezi/-pdf-extractor-app-.git
cd pdf-extractor-app
```

2. Install requirements:
```bash
pip install -r requirements.txt
```

3. Run the app:
```bash
streamlit run src/agent/pdf_extractor_streamlit.py
```

## Configuration

1. Get an OpenRouter API key from https://openrouter.ai/keys
2. Enter the API key in the app's sidebar

## Project Structure

```
pdf-extractor-app/
├── src/
│   └── agent/
│       ├── pdf_extractor_streamlit.py
│       ├── ai_processor.py
│       └── math_example.py
├── requirements.txt
└── README.md
```

## Requirements

- Python 3.8+
- Streamlit
- Pandas
- PyPDF2
- Openpyxl
- Requests

## License

MIT License
