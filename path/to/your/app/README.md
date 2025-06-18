# PDF Extractor Application

A powerful application for extracting data from PDF files, available in both GUI and Streamlit versions.

## Features

- Extract data from multiple PDF files simultaneously
- Two interface options:
  - Desktop GUI application
  - Web-based Streamlit application
- Save extracted data to Excel files
- Customizable output directory
- Progress tracking and error handling
- Modern and user-friendly interface

## Installation

1. Clone the repository:
```bash
git clone https://github.com/gupadovezi/-pdf-extractor-app-.git
cd -pdf-extractor-app-
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

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
│       └── math_example.py           # PDF processing logic
├── requirements.txt                  # Project dependencies
└── README.md                        # This file
```

## Requirements

- Python 3.8+
- Dependencies listed in requirements.txt:
  - streamlit==1.45.1
  - pandas==2.3.0
  - PyPDF2==3.0.1
  - openpyxl==3.1.5

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
