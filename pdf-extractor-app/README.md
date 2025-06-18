# AI-Powered PDF Extractor Application

A powerful application for extracting and analyzing data from PDF files using AI, available in both GUI and Streamlit versions.

## Features

- Extract data from multiple PDF files simultaneously
- AI-powered analysis of research papers
- Two interface options:
  - Desktop GUI application
  - Web-based Streamlit application
- Save extracted data to Excel files
- Customizable output directory
- Progress tracking and error handling
- Modern and user-friendly interface
- AI-powered insights and analysis

## AI Capabilities

The application uses OpenAI's GPT models to:
- Extract structured information from PDFs
- Analyze research findings
- Identify patterns and insights
- Generate comprehensive reports
- Compare multiple papers
- Identify research gaps and future directions

## Installation

1. Clone the repository:
```bash
git clone https://github.com/gupadovezi/-pdf-extractor-app-.git
cd pdf-extractor-app
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Set up your OpenAI API key:
   - Get an API key from [OpenAI](https://platform.openai.com/api-keys)
   - Create a `.env` file in the project root
   - Add your API key: `OPENAI_API_KEY=your_api_key_here`

## Usage

### GUI Version
Run the desktop application:
```bash
python src/agent/pdf_extractor_gui.py
```

### Streamlit Version
Run the web application:
```bash
streamlit run src/agent/pdf_extractor_streamlit.py
```

## Project Structure

```
pdf-extractor-app/
├── src/
│   └── agent/
│       ├── pdf_extractor_gui.py      # Desktop GUI application
│       ├── pdf_extractor_streamlit.py # Web application
│       ├── math_example.py           # PDF processing logic
│       └── ai_processor.py           # AI processing functionality
├── requirements.txt                  # Project dependencies
└── README.md                        # This file
```

## Requirements

- Python 3.8+
- OpenAI API key
- Dependencies listed in requirements.txt:
  - streamlit==1.45.1
  - pandas==2.3.0
  - PyPDF2==3.0.1
  - openpyxl==3.1.5
  - openai==1.12.0
  - python-dotenv==1.0.0
  - langchain==0.1.0
  - langchain-openai==0.0.2

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
