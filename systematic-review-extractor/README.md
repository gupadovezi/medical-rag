# Systematic Review Extractor with RAG

This application helps you extract and analyze information from PDF files using AI and RAG (Retrieval-Augmented Generation). It uses OpenRouter's API to process research papers, build a searchable knowledge base, and provide contextual answers to questions.

## Features

- Upload multiple PDF files for processing
- Build a RAG knowledge base from uploaded PDFs
- Extract key information using AI
- Generate Excel file with structured data
- Analyze findings across multiple papers
- Ask questions and get contextual answers from your research papers
- Identify common themes and research gaps
- Vector-based document search and retrieval

## Online Version

The app is deployed at: [https://wwwsystematicreviewextractor.streamlit.app](https://wwwsystematicreviewextractor.streamlit.app)

## How RAG Works

1. **Document Processing**: PDFs are processed and split into chunks
2. **Embedding Generation**: Each chunk is converted to vector embeddings
3. **Knowledge Base**: Documents are stored in a searchable vector database
4. **Query Processing**: Questions are converted to embeddings and matched with relevant documents
5. **Contextual Answers**: AI generates answers using retrieved context

## Deployment

### Deploying to Streamlit Cloud

1. Fork this repository to your GitHub account
2. Go to [share.streamlit.io](https://share.streamlit.io/)
3. Sign in with your GitHub account
4. Click "New app"
5. Select your forked repository
6. Configure the app:
   - Main file path: `app.py`
   - Python version: 3.8 or higher
7. Add your OpenRouter API key in the "Advanced settings":
   - Click "Add secret"
   - Key: `OPENROUTER_API_KEY`
   - Value: Your OpenRouter API key
8. Click "Deploy"

### Local Development

1. Clone the repository:
```bash
git clone https://github.com/yourusername/systematic-review-extractor.git
cd systematic-review-extractor
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
streamlit run app.py
```

## Usage

### Building the Knowledge Base

1. Upload PDF files using the sidebar "Upload PDFs for Knowledge Base"
2. Click "Add to Knowledge Base" to process and embed the documents
3. Monitor the knowledge base stats to see progress

### Asking Questions

1. Use the "RAG Query" section in the sidebar
2. Type your question about the research papers
3. Click "Ask RAG" to get contextual answers
4. View relevant documents and AI-generated responses

### Processing PDFs

1. Upload PDF files in the main area
2. Click "Process PDFs" to extract structured information
3. Download Excel files with extracted data

## Configuration

1. Get an OpenRouter API key from https://openrouter.ai/keys
2. For local development: Add the API key to `.streamlit/secrets.toml`
3. For Streamlit Cloud: Add the API key in the app's secrets section

## Project Structure

```
systematic-review-extractor/
├── .streamlit/
│   └── secrets.toml
├── app.py
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
- NumPy

## License

MIT License 