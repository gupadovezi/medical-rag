# RAG-Enhanced Chain Implementation

This module demonstrates how to enhance the basic LLM chain with RAG (Retrieval-Augmented Generation) capabilities using LangGraph and LangChain.

## What is RAG?

RAG (Retrieval-Augmented Generation) is a technique that combines:
- **Retrieval**: Finding relevant documents from a knowledge base
- **Augmentation**: Adding retrieved context to the LLM's input
- **Generation**: Using the enhanced context to generate more accurate responses

## Files in this Module

### 1. `rag-chain.ipynb`
A comprehensive Jupyter notebook that demonstrates:
- Document storage with vector embeddings
- Retrieval system setup
- RAG tools creation
- Enhanced state management
- Interactive testing

### 2. `rag_chain.py`
A complete Python script that includes:
- Full RAG system implementation
- Interactive chat interface
- Comprehensive testing suite
- Error handling and user experience features

### 3. `rag_demo.py`
A quick demo script for:
- Fast testing of RAG functionality
- Simple example with minimal setup
- Understanding the core concepts

## Quick Start

### Prerequisites

1. **Install Dependencies**
   ```bash
   pip install langchain_openai langchain_core langgraph langchain_community chromadb sentence_transformers
   ```

2. **Set Environment Variables**
   ```bash
   export OPENAI_API_KEY="your-openai-api-key-here"
   ```

### Running the Demo

1. **Quick Demo** (Recommended for first-time users):
   ```bash
   python module-1/rag_demo.py
   ```

2. **Full RAG System**:
   ```bash
   python module-1/rag_chain.py
   ```

3. **Jupyter Notebook**:
   ```bash
   jupyter notebook module-1/rag-chain.ipynb
   ```

## Key Features

### 🔍 Document Retrieval
- Vector-based similarity search using OpenAI embeddings
- Configurable number of retrieved documents
- Metadata tracking for source attribution

### 🛠️ RAG Tools
- `search_knowledge_base`: Search for relevant information
- `list_available_topics`: List all available topics
- `summarize_documents`: Generate summaries with sources

### 🧠 Enhanced State Management
- Extended `MessagesState` with context tracking
- Automatic context retrieval for user queries
- Tool call integration

### 🔄 Graph-Based Architecture
- LangGraph integration for complex workflows
- Conditional document retrieval
- Scalable architecture for production use

## Architecture Overview

```
User Query → Document Retrieval → Context Augmentation → LLM Generation → Response
     ↓              ↓                      ↓                    ↓
  Vector Store → Similarity Search → Context Addition → Enhanced Response
```

## Example Usage

### Basic RAG Query
```python
from langchain_core.messages import HumanMessage

# Query the RAG system
result = rag_graph.invoke({
    "messages": [HumanMessage(content="What is machine learning?")]
})

# Get the response
response = result['messages'][-1].content
print(response)
```

### Using RAG Tools
```python
# Search the knowledge base
search_result = search_knowledge_base("deep learning frameworks")

# List available topics
topics = list_available_topics()

# Summarize documents
summary = summarize_documents("AI applications in healthcare")
```

## Customization

### Adding Your Own Documents
```python
from langchain_core.documents import Document

# Create your documents
my_documents = [
    Document(
        page_content="Your document content here...",
        metadata={"source": "my_doc.txt", "topic": "my_topic"}
    )
]

# Add to vector store
vectorstore = Chroma.from_documents(
    documents=my_documents,
    embedding=embeddings,
    collection_name="my_knowledge_base"
)
```

### Custom Retrieval Strategy
```python
def custom_retrieval(query: str, k: int = 5):
    """Custom retrieval with different parameters."""
    docs = vectorstore.similarity_search(
        query, 
        k=k,
        filter={"topic": "specific_topic"}  # Add filters
    )
    return docs
```

### Enhanced RAG Tools
```python
def advanced_search(query: str, filters: dict = None) -> str:
    """Advanced search with filtering capabilities."""
    # Implementation here
    pass

# Add to LLM tools
llm_with_tools = llm.bind_tools([advanced_search, other_tools])
```

## Best Practices

### 1. Document Quality
- Use high-quality, relevant documents
- Include proper metadata for filtering
- Chunk documents appropriately for retrieval

### 2. Retrieval Optimization
- Tune the number of retrieved documents (k parameter)
- Use appropriate similarity metrics
- Consider hybrid search strategies

### 3. Context Management
- Limit context length to avoid token limits
- Prioritize most relevant information
- Include source attribution

### 4. Error Handling
- Handle cases where no relevant documents are found
- Implement fallback strategies
- Provide meaningful error messages

## Troubleshooting

### Common Issues

1. **API Key Not Set**
   ```
   ❌ Please set your OPENAI_API_KEY environment variable
   ```
   **Solution**: Set the environment variable or use the interactive prompt

2. **Vector Store Errors**
   ```
   ❌ Error creating vector store
   ```
   **Solution**: Check internet connection and API key validity

3. **Memory Issues**
   ```
   ❌ Out of memory
   ```
   **Solution**: Reduce document size or use smaller embedding models

### Performance Tips

- Use smaller embedding models for faster retrieval
- Implement caching for frequently accessed documents
- Consider using async operations for large document sets
- Monitor token usage to optimize costs

## Advanced Features

### Multi-Modal RAG
Extend the system to handle images, audio, and other media types.

### Hybrid Search
Combine vector search with keyword-based search for better results.

### Conversation Memory
Implement long-term memory for multi-turn conversations.

### Real-time Updates
Add capability to update the knowledge base in real-time.

## Contributing

To extend this RAG system:

1. Fork the repository
2. Create a feature branch
3. Implement your enhancements
4. Add tests and documentation
5. Submit a pull request

## Resources

- [LangChain Documentation](https://python.langchain.com/)
- [LangGraph Documentation](https://langchain-ai.github.io/langgraph/)
- [Chroma Vector Store](https://docs.trychroma.com/)
- [OpenAI Embeddings](https://platform.openai.com/docs/guides/embeddings)

## License

This implementation is part of the LangChain Academy and follows the same licensing terms. 