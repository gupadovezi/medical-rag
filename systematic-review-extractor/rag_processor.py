import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import json
import tempfile
import os
from datetime import datetime
import requests
from typing import List, Dict, Any
import PyPDF2
import io

class RAGProcessor:
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.embeddings_cache = {}
        self.knowledge_base = []
        
    def extract_text_from_pdf(self, pdf_file) -> str:
        """Extract text from a PDF file."""
        try:
            pdf_reader = PyPDF2.PdfReader(pdf_file)
            text = ""
            for page in pdf_reader.pages:
                text += page.extract_text() + "\n"
            return text
        except Exception as e:
            st.error(f"Error extracting text from PDF: {str(e)}")
            return ""
    
    def chunk_text(self, text: str, chunk_size: int = 1000, overlap: int = 200) -> List[str]:
        """Split text into overlapping chunks."""
        chunks = []
        start = 0
        while start < len(text):
            end = start + chunk_size
            chunk = text[start:end]
            chunks.append(chunk)
            start = end - overlap
        return chunks
    
    def get_embeddings(self, text: str) -> List[float]:
        """Get embeddings for text using OpenRouter API."""
        url = "https://openrouter.ai/api/v1/embeddings"
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "HTTP-Referer": "https://wwwsystematicreviewextractor.streamlit.app",
            "Content-Type": "application/json"
        }
        payload = {
            "model": "text-embedding-ada-002",
            "input": text
        }
        
        try:
            response = requests.post(url, headers=headers, json=payload)
            response.raise_for_status()
            return response.json()['data'][0]['embedding']
        except Exception as e:
            st.error(f"Error getting embeddings: {str(e)}")
            return []
    
    def add_to_knowledge_base(self, pdf_files: List, metadata: Dict = None) -> None:
        """Add PDF files to the knowledge base."""
        with st.spinner("Processing PDFs for knowledge base..."):
            for pdf_file in pdf_files:
                # Extract text
                text = self.extract_text_from_pdf(pdf_file)
                if not text:
                    continue
                
                # Chunk text
                chunks = self.chunk_text(text)
                
                # Get embeddings for each chunk
                for i, chunk in enumerate(chunks):
                    embedding = self.get_embeddings(chunk)
                    if embedding:
                        self.knowledge_base.append({
                            'text': chunk,
                            'embedding': embedding,
                            'metadata': {
                                'filename': pdf_file.name,
                                'chunk_id': i,
                                'timestamp': datetime.now().isoformat(),
                                **(metadata or {})
                            }
                        })
        
        st.success(f"Added {len(pdf_files)} files to knowledge base")
    
    def search_knowledge_base(self, query: str, top_k: int = 5) -> List[Dict]:
        """Search the knowledge base for relevant documents."""
        if not self.knowledge_base:
            return []
        
        # Get query embedding
        query_embedding = self.get_embeddings(query)
        if not query_embedding:
            return []
        
        # Calculate similarities
        similarities = []
        for doc in self.knowledge_base:
            similarity = self.cosine_similarity(query_embedding, doc['embedding'])
            similarities.append((similarity, doc))
        
        # Sort by similarity and return top_k
        similarities.sort(key=lambda x: x[0], reverse=True)
        return [doc for _, doc in similarities[:top_k]]
    
    def cosine_similarity(self, vec1: List[float], vec2: List[float]) -> float:
        """Calculate cosine similarity between two vectors."""
        vec1 = np.array(vec1)
        vec2 = np.array(vec2)
        return np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    
    def rag_query(self, query: str, context_window: int = 3000) -> str:
        """Perform RAG query with context from knowledge base."""
        # Search for relevant documents
        relevant_docs = self.search_knowledge_base(query)
        
        if not relevant_docs:
            return "No relevant documents found in knowledge base."
        
        # Build context from relevant documents
        context = ""
        for doc in relevant_docs:
            context += f"Document: {doc['metadata']['filename']}\n"
            context += f"Content: {doc['text']}\n\n"
            if len(context) > context_window:
                break
        
        # Create enhanced prompt with context
        enhanced_prompt = f"""
        Based on the following context from research papers, answer the question.
        
        Context:
        {context}
        
        Question: {query}
        
        Please provide a comprehensive answer based on the context provided. If the context doesn't contain enough information to answer the question, please state that clearly.
        """
        
        return enhanced_prompt
    
    def save_knowledge_base(self, filepath: str) -> None:
        """Save knowledge base to file."""
        try:
            # Convert numpy arrays to lists for JSON serialization
            serializable_kb = []
            for doc in self.knowledge_base:
                serializable_doc = {
                    'text': doc['text'],
                    'embedding': doc['embedding'],
                    'metadata': doc['metadata']
                }
                serializable_kb.append(serializable_doc)
            
            with open(filepath, 'w') as f:
                json.dump(serializable_kb, f, indent=2)
            st.success(f"Knowledge base saved to {filepath}")
        except Exception as e:
            st.error(f"Error saving knowledge base: {str(e)}")
    
    def load_knowledge_base(self, filepath: str) -> None:
        """Load knowledge base from file."""
        try:
            with open(filepath, 'r') as f:
                self.knowledge_base = json.load(f)
            st.success(f"Knowledge base loaded from {filepath}")
        except Exception as e:
            st.error(f"Error loading knowledge base: {str(e)}")
    
    def get_knowledge_base_stats(self) -> Dict[str, Any]:
        """Get statistics about the knowledge base."""
        if not self.knowledge_base:
            return {"total_documents": 0, "total_chunks": 0, "unique_files": 0}
        
        unique_files = set(doc['metadata']['filename'] for doc in self.knowledge_base)
        return {
            "total_documents": len(unique_files),
            "total_chunks": len(self.knowledge_base),
            "unique_files": len(unique_files)
        } 