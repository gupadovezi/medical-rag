#!/usr/bin/env python3
"""
RAG-Enhanced Chain Implementation

This script demonstrates how to enhance a basic LLM chain with RAG (Retrieval-Augmented Generation) capabilities.
It includes document storage, retrieval, and integration with LangGraph chains.
"""

import os
import getpass
from typing import List, Dict, Any
from langchain_core.messages import HumanMessage, AIMessage, SystemMessage
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain_community.vectorstores import Chroma
from langchain_core.documents import Document
from langgraph.graph import StateGraph, MessagesState, START, END
from typing_extensions import Annotated
from langgraph.graph.message import add_messages


def setup_environment():
    """Set up environment variables for API keys."""
    def _set_env(var: str):
        if not os.environ.get(var):
            os.environ[var] = getpass.getpass(f"{var}: ")
    
    _set_env("OPENAI_API_KEY")
    print("✅ Environment setup complete")


def create_sample_documents():
    """Create sample documents for demonstration."""
    sample_documents = [
        Document(
            page_content="""
            Artificial Intelligence (AI) is a branch of computer science that aims to create intelligent machines 
            that can perform tasks that typically require human intelligence. These tasks include learning, 
            reasoning, problem-solving, perception, and language understanding. AI has applications in various 
            fields including healthcare, finance, transportation, and entertainment.
            """,
            metadata={"source": "ai_intro.txt", "topic": "artificial_intelligence"}
        ),
        Document(
            page_content="""
            Machine Learning is a subset of AI that enables computers to learn and improve from experience 
            without being explicitly programmed. It uses algorithms to identify patterns in data and make 
            predictions or decisions. Common types include supervised learning, unsupervised learning, and 
            reinforcement learning.
            """,
            metadata={"source": "ml_basics.txt", "topic": "machine_learning"}
        ),
        Document(
            page_content="""
            Deep Learning is a subset of machine learning that uses neural networks with multiple layers 
            to model and understand complex patterns. It has revolutionized fields like computer vision, 
            natural language processing, and speech recognition. Popular frameworks include TensorFlow and PyTorch.
            """,
            metadata={"source": "deep_learning.txt", "topic": "deep_learning"}
        ),
        Document(
            page_content="""
            Natural Language Processing (NLP) is a field of AI that focuses on the interaction between 
            computers and human language. It enables machines to understand, interpret, and generate human 
            language. Applications include chatbots, translation services, and sentiment analysis.
            """,
            metadata={"source": "nlp_overview.txt", "topic": "natural_language_processing"}
        ),
        Document(
            page_content="""
            Large Language Models (LLMs) are a type of neural network that can understand and generate 
            human language. They are trained on vast amounts of text data and can perform tasks like 
            text generation, translation, summarization, and question answering. Examples include GPT, 
            BERT, and T5 models.
            """,
            metadata={"source": "llm_guide.txt", "topic": "large_language_models"}
        )
    ]
    
    print(f"✅ Created {len(sample_documents)} sample documents")
    return sample_documents


def setup_vector_store(documents: List[Document]):
    """Set up vector store with documents."""
    print("🔧 Setting up vector store...")
    
    # Initialize embeddings
    embeddings = OpenAIEmbeddings()
    
    # Create vector store
    vectorstore = Chroma.from_documents(
        documents=documents,
        embedding=embeddings,
        collection_name="ai_knowledge_base"
    )
    
    print(f"✅ Vector store created with {vectorstore._collection.count()} documents")
    return vectorstore


def retrieve_documents(vectorstore: Chroma, query: str, k: int = 3) -> List[Document]:
    """Retrieve relevant documents from the vector store."""
    docs = vectorstore.similarity_search(query, k=k)
    return docs


def search_knowledge_base(vectorstore: Chroma, query: str) -> str:
    """Search the knowledge base for relevant information."""
    docs = retrieve_documents(vectorstore, query, k=3)
    
    if not docs:
        return "No relevant information found in the knowledge base."
    
    result = "Relevant information from knowledge base:\n\n"
    for i, doc in enumerate(docs, 1):
        result += f"Document {i} (Source: {doc.metadata['source']}):\n"
        result += f"{doc.page_content.strip()}\n\n"
    
    return result


def list_available_topics(documents: List[Document]) -> str:
    """List all available topics in the knowledge base."""
    topics = [doc.metadata['topic'] for doc in documents]
    return f"Available topics in knowledge base: {', '.join(topics)}"


def summarize_documents(vectorstore: Chroma, query: str) -> str:
    """Retrieve and summarize relevant documents for a query."""
    docs = retrieve_documents(vectorstore, query, k=3)
    
    if not docs:
        return "No relevant information found."
    
    # Create a summary prompt
    summary_prompt = f"""
    Based on the following documents, provide a comprehensive summary for the query: '{query}'
    
    Documents:
    """
    
    for i, doc in enumerate(docs, 1):
        summary_prompt += f"\nDocument {i} (Source: {doc.metadata['source']}):\n{doc.page_content}\n"
    
    summary_prompt += "\nPlease provide a comprehensive summary with proper attribution to sources."
    
    # Use the LLM to generate a summary
    summary_llm = ChatOpenAI(model="gpt-4o")
    summary = summary_llm.invoke([HumanMessage(content=summary_prompt)])
    
    return summary.content


class RAGState(MessagesState):
    """Extended state that includes context from retrieved documents."""
    context: Annotated[str, add_messages] = ""


def create_rag_tools(vectorstore: Chroma, documents: List[Document]):
    """Create RAG tools with access to the vector store."""
    def search_tool(query: str) -> str:
        return search_knowledge_base(vectorstore, query)
    
    def list_topics_tool() -> str:
        return list_available_topics(documents)
    
    def summarize_tool(query: str) -> str:
        return summarize_documents(vectorstore, query)
    
    return [search_tool, list_topics_tool, summarize_tool]


def create_rag_graph(vectorstore: Chroma):
    """Create the RAG-enhanced graph."""
    print("🔧 Creating RAG-enhanced graph...")
    
    # Initialize LLM with RAG tools
    llm = ChatOpenAI(model="gpt-4o")
    rag_tools = create_rag_tools(vectorstore, create_sample_documents())
    llm_with_rag_tools = llm.bind_tools(rag_tools)
    
    def rag_llm_node(state: RAGState):
        """Node that processes messages with RAG capabilities."""
        # Get the last message
        last_message = state["messages"][-1]
        
        # If it's a human message, try to retrieve relevant context
        if isinstance(last_message, HumanMessage):
            # Retrieve relevant documents
            docs = retrieve_documents(vectorstore, last_message.content, k=2)
            
            # Create context from retrieved documents
            context = ""
            if docs:
                context = "Relevant context:\n"
                for i, doc in enumerate(docs, 1):
                    context += f"\nDocument {i} (Source: {doc.metadata['source']}):\n"
                    context += f"{doc.page_content.strip()}\n"
            
            # Add context to state
            if context:
                state["context"] = context
        
        # Invoke LLM with tools
        result = llm_with_rag_tools.invoke(state["messages"])
        
        return {"messages": [result]}
    
    # Build the RAG-enhanced graph
    builder = StateGraph(RAGState)
    builder.add_node("rag_llm", rag_llm_node)
    builder.add_edge(START, "rag_llm")
    builder.add_edge("rag_llm", END)
    rag_graph = builder.compile()
    
    print("✅ RAG-enhanced graph created successfully")
    return rag_graph


def test_rag_system(rag_graph):
    """Test the RAG system with various queries."""
    print("\n🧪 Testing RAG System")
    print("=" * 50)
    
    test_queries = [
        "Hello! How are you?",
        "What is artificial intelligence?",
        "Explain machine learning and its types",
        "What topics are available in your knowledge base?",
        "Tell me about deep learning frameworks"
    ]
    
    for i, query in enumerate(test_queries, 1):
        print(f"\nTest {i}: {query}")
        print("-" * 30)
        
        messages = rag_graph.invoke({"messages": [HumanMessage(content=query)]})
        
        # Get the AI response
        ai_message = messages['messages'][-1]
        print(f"🤖 Assistant: {ai_message.content}")
        
        # Show tool calls if any
        if hasattr(ai_message, 'tool_calls') and ai_message.tool_calls:
            print("🔧 Tool calls made:")
            for tool_call in ai_message.tool_calls:
                print(f"  - {tool_call['name']}")
        
        print()


def interactive_chat(rag_graph):
    """Interactive chat function."""
    print("\n🤖 RAG-Enhanced AI Assistant")
    print("=" * 50)
    print("I have access to a knowledge base about AI, ML, and related technologies.")
    print("You can ask me questions about these topics, and I'll retrieve relevant information.")
    print("Type 'quit' to exit.\n")
    
    # Initialize conversation state
    conversation_state = RAGState(messages=[], context="")
    
    while True:
        try:
            user_input = input("\n👤 You: ")
            
            if user_input.lower() in ['quit', 'exit', 'bye']:
                print("\n🤖 Assistant: Goodbye! Have a great day!")
                break
            
            # Add user message to state
            conversation_state["messages"].append(HumanMessage(content=user_input))
            
            # Process with RAG graph
            result = rag_graph.invoke(conversation_state)
            
            # Get the AI response
            ai_message = result["messages"][-1]
            
            # Update conversation state
            conversation_state = result
            
            # Display response
            print(f"\n🤖 Assistant: {ai_message.content}")
            
            # If there are tool calls, show them
            if hasattr(ai_message, 'tool_calls') and ai_message.tool_calls:
                print("\n🔧 Tool calls made:")
                for tool_call in ai_message.tool_calls:
                    print(f"  - {tool_call['name']}")
                    
        except KeyboardInterrupt:
            print("\n\n🤖 Assistant: Goodbye! Have a great day!")
            break
        except Exception as e:
            print(f"\n❌ Error: {e}")
            print("Please try again.")


def main():
    """Main function to run the RAG system."""
    print("🚀 Starting RAG-Enhanced Chain System")
    print("=" * 50)
    
    # Setup environment
    setup_environment()
    
    # Create sample documents
    documents = create_sample_documents()
    
    # Setup vector store
    vectorstore = setup_vector_store(documents)
    
    # Create RAG graph
    rag_graph = create_rag_graph(vectorstore)
    
    # Test the system
    test_rag_system(rag_graph)
    
    # Start interactive chat
    print("\n" + "=" * 50)
    print("Starting interactive chat mode...")
    interactive_chat(rag_graph)


if __name__ == "__main__":
    main() 