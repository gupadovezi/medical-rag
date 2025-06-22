![LangChain Academy](https://cdn.prod.website-files.com/65b8cd72835ceeacd4449a53/66e9eba1020525eea7873f96_LCA-big-green%20(2).svg)

## Introduction

Welcome to LangChain Academy! 
This is a growing set of modules focused on foundational concepts within the LangChain ecosystem. 
Module 0 is basic setup and Modules 1 - 4 focus on LangGraph, progressively adding more advanced themes. 
In each module folder, you'll see a set of notebooks. A LangChain Academy accompanies each notebook 
to guide you through the topic. Each module also has a `studio` subdirectory, with a set of relevant 
graphs that we will explore using the LangGraph API and Studio.

## Setup

### Python version

To get the most out of this course, please ensure you're using Python 3.11 or later. 
This version is required for optimal compatibility with LangGraph. If you're on an older version, 
upgrading will ensure everything runs smoothly.
```
python3 --version
```

### Clone repo
```
git clone https://github.com/langchain-ai/langchain-academy.git
$ cd langchain-academy
```

### Create an environment and install dependencies
#### Mac/Linux/WSL
```
$ python3 -m venv lc-academy-env
$ source lc-academy-env/bin/activate
$ pip install -r requirements.txt
```
#### Windows Powershell
```
PS> python3 -m venv lc-academy-env
PS> Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process
PS> lc-academy-env\scripts\activate
PS> pip install -r requirements.txt
```

### Running notebooks
If you don't have Jupyter set up, follow installation instructions [here](https://jupyter.org/install).
```
$ jupyter notebook
```

### Setting up env variables
Briefly going over how to set up environment variables. You can also 
use a `.env` file with `python-dotenv` library.
#### Mac/Linux/WSL
```
$ export API_ENV_VAR="your-api-key-here"
```
#### Windows Powershell
```
PS> $env:API_ENV_VAR = "your-api-key-here"
```

### Set OpenAI API key
* If you don't have an OpenAI API key, you can sign up [here](https://openai.com/index/openai-api/).
*  Set `OPENAI_API_KEY` in your environment 

### Sign up and Set LangSmith API
* Sign up for LangSmith [here](https://smith.langchain.com/), find out more about LangSmith
* and how to use it within your workflow [here](https://www.langchain.com/langsmith), and relevant library [docs](https://docs.smith.langchain.com/)!
*  Set `LANGCHAIN_API_KEY`, `LANGCHAIN_TRACING_V2=true` in your environment 

### Set up Tavily API for web search

* Tavily Search API is a search engine optimized for LLMs and RAG, aimed at efficient, 
quick, and persistent search results. 
* You can sign up for an API key [here](https://tavily.com/). 
It's easy to sign up and offers a very generous free tier. Some lessons (in Module 4) will use Tavily. 

* Set `TAVILY_API_KEY` in your environment.

### Set up LangGraph Studio

* LangGraph Studio is a custom IDE for viewing and testing agents.
* Studio can be run locally and opened in your browser on Mac, Windows, and Linux.
* See documentation [here](https://langchain-ai.github.io/langgraph/concepts/langgraph_studio/#local-development-server) on the local Studio development server and [here](https://langchain-ai.github.io/langgraph/how-tos/local-studio/#run-the-development-server). 
* Graphs for LangGraph Studio are in the `module-x/studio/` folders.
* To start the local development server, run the following command in your terminal in the `/studio` directory each module:

```
langgraph dev
```

You should see the following output:
```
- 🚀 API: http://127.0.0.1:2024
- 🎨 Studio UI: https://smith.langchain.com/studio/?baseUrl=http://127.0.0.1:2024
- 📚 API Docs: http://127.0.0.1:2024/docs
```

Open your browser and navigate to the Studio UI: `https://smith.langchain.com/studio/?baseUrl=http://127.0.0.1:2024`.

* To use Studio, you will need to create a .env file with the relevant API keys
* Run this from the command line to create these files for module 1 to 6, as an example:
```
for i in {1..6}; do
  cp module-$i/studio/.env.example module-$i/studio/.env
  echo "OPENAI_API_KEY=\"$OPENAI_API_KEY\"" > module-$i/studio/.env
done
echo "TAVILY_API_KEY=\"$TAVILY_API_KEY\"" >> module-4/studio/.env
```

# 🏥 Medical RAG Assistant

A comprehensive medical Retrieval-Augmented Generation (RAG) system with a modern web interface for musculoskeletal pain conditions.

## 🚀 Live Demo

**Web App**: [Deploy to Streamlit Cloud](https://share.streamlit.io)

## ✨ Features

- **🤖 AI-Powered Medical Assistant** - RAG system with medical knowledge base
- **📚 Comprehensive Knowledge Base** - Covers 10+ musculoskeletal conditions
- **🔬 PubMed Integration** - Real-time medical research data
- **💬 Interactive Chat Interface** - Modern Streamlit web UI
- **🔄 Offline Mode** - Works without API calls
- **📊 Chat History** - Persistent conversation tracking
- **⚡ Rate Limit Handling** - Automatic fallback and retry logic

## 🏗️ Architecture

- **Frontend**: Streamlit web interface
- **Backend**: Python with LangChain
- **Vector Database**: ChromaDB
- **Embeddings**: HuggingFace Transformers
- **Medical Data**: PubMed API + Custom knowledge base

## 🛠️ Quick Start

### Local Development

1. **Clone the repository**:
   ```bash
   git clone https://github.com/gupadovezi/-pdf-extractor-app-.git
   cd langchain-academy-main
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the web app**:
   ```bash
   streamlit run rag_web_app.py
   ```

4. **Open your browser**: http://localhost:8501

### Deployment

See [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) for detailed deployment instructions.

## 📁 Project Structure

```
├── rag_web_app.py              # Main Streamlit web application
├── run_rag_interactive.py      # Terminal-based interactive version
├── run_rag_offline.py          # Offline-only version
├── requirements.txt            # Python dependencies
├── .streamlit/config.toml      # Streamlit configuration
└── DEPLOYMENT_GUIDE.md         # Deployment instructions
```

## 🎯 Use Cases

- **Medical Research**: Query medical literature and research papers
- **Patient Education**: Get information about musculoskeletal conditions
- **Clinical Support**: Access evidence-based medical information
- **Medical Training**: Educational tool for healthcare professionals

## 🔧 Configuration

The app supports multiple modes:
- **Online Mode**: Uses OpenRouter API for enhanced responses
- **Offline Mode**: Works entirely with local knowledge base
- **Hybrid Mode**: Combines local knowledge with API when available

## 📊 Performance

- **Response Time**: < 5 seconds for most queries
- **Knowledge Base**: 10+ musculoskeletal conditions covered
- **PubMed Integration**: Real-time access to medical research
- **Offline Capability**: 100% functional without internet

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## 📄 License

This project is licensed under the MIT License.

## 🆘 Support

For issues and questions:
- Create a GitHub issue
- Check the deployment guide
- Review the documentation

---

**Built with ❤️ using Streamlit, LangChain, and medical research data**
