#!/usr/bin/env python3
"""
Quick RAG Demo

This script provides a quick demonstration of the RAG-enhanced chain system.
"""

import os
from langchain_core.messages import HumanMessage
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain_community.vectorstores import Chroma
from langchain_core.documents import Document
from langgraph.graph import StateGraph, MessagesState, START, END
from typing_extensions import Annotated
from langgraph.graph.message import add_messages


def quick_rag_demo():
    """Quick demonstration of RAG functionality."""
    
    # Check for API key
    if not os.environ.get("OPENAI_API_KEY"):
        print("❌ Please set your OPENAI_API_KEY environment variable")
        return
    
    print("🚀 Quick RAG Demo")
    print("=" * 40)
    
    # Create a simple document
    documents = [
        Document(
            page_content="""
            Python is a high-level, interpreted programming language known for its simplicity and readability.
            It was created by Guido van Rossum and first released in 1991. Python supports multiple programming
            paradigms including procedural, object-oriented, and functional programming. It's widely used in
            web development, data science, artificial intelligence, and automation.
            """,
            metadata={"source": "python_intro.txt", "topic": "programming"}
        )
    ]
    
    # Setup vector store
    print("🔧 Setting up vector store...")
    embeddings = OpenAIEmbeddings()
    vectorstore = Chroma.from_documents(
        documents=documents,
        embedding=embeddings,
        collection_name="demo_knowledge_base"
    )
    
    # Create RAG state
    class RAGState(MessagesState):
        context: Annotated[str, add_messages] = ""
    
    # Create LLM with RAG capabilities
    llm = ChatOpenAI(model="gpt-4o")
    
    def rag_node(state: RAGState):
        """Simple RAG node that retrieves context and generates response."""
        last_message = state["messages"][-1]
        
        if isinstance(last_message, HumanMessage):
            # Retrieve relevant documents
            docs = vectorstore.similarity_search(last_message.content, k=1)
            
            # Create context
            context = ""
            if docs:
                context = f"Relevant information: {docs[0].page_content.strip()}"
                state["context"] = context
        
        # Generate response
        result = llm.invoke(state["messages"])
        return {"messages": [result]}
    
    # Build graph
    builder = StateGraph(RAGState)
    builder.add_node("rag_llm", rag_node)
    builder.add_edge(START, "rag_llm")
    builder.add_edge("rag_llm", END)
    graph = builder.compile()
    
    print("✅ RAG system ready!")
    
    # Test queries
    test_queries = [
        "What is Python?",
        "Who created Python?",
        "What are Python's main features?",
        "Hello, how are you?"
    ]
    
    print("\n🧪 Testing RAG System:")
    print("-" * 30)
    
    for i, query in enumerate(test_queries, 1):
        print(f"\nTest {i}: {query}")
        print("-" * 20)
        
        result = graph.invoke({"messages": [HumanMessage(content=query)]})
        ai_message = result['messages'][-1]
        
        print(f"🤖 Assistant: {ai_message.content}")
        
        # Show context if available
        if result.get('context'):
            print(f"📚 Context used: {result['context'][:100]}...")
    
    print("\n✅ Demo completed!")


if __name__ == "__main__":
    quick_rag_demo() 