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

## Deployment

### Deploying to Streamlit Cloud

1. Fork this repository to your GitHub account
2. Go to [share.streamlit.io](https://share.streamlit.io/)
3. Sign in with your GitHub account
4. Click "New app"
5. Select your forked repository
6. Configure the app:
   - Main file path: `src/agent/pdf_extractor_streamlit.py`
   - Python version: 3.8 or higher
7. Add your OpenRouter API key in the "Advanced settings":
   - Click "Add secret"
   - Key: `OPENROUTER_API_KEY`
   - Value: Your OpenRouter API key
8. Click "Deploy"

### Local Development

1. Clone the repository:
```bash
git clone https://github.com/gupadovezi/-pdf-extractor-app-.git
cd pdf-extractor-app
```

2. Install requirements:
```bash
pip install -r requirements.txt
```

3. Create a `.streamlit/secrets.toml` file:
```toml
OPENROUTER_API_KEY = "your-api-key-here"
```

4. Run the app:
```bash
streamlit run src/agent/pdf_extractor_streamlit.py
```

## Configuration

1. Get an OpenRouter API key from https://openrouter.ai/keys
2. For local development: Add the API key to `.streamlit/secrets.toml`
3. For Streamlit Cloud: Add the API key in the app's secrets section

## Project Structure

```
pdf-extractor-app/
├── .streamlit/
│   └── secrets.toml
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
