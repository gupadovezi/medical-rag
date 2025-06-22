import streamlit as st
import pandas as pd
from pathlib import Path
import os
from datetime import datetime
import tempfile
import shutil
import requests
import json
import numpy as np
from typing import List, Dict, Any
import PyPDF2

# Constants
MAX_FILE_SIZE = 10 * 1024 * 1024  # 10MB
ALLOWED_EXTENSIONS = {'pdf'}

# Set page config
st.set_page_config(
    page_title="Systematic Review Extractor with RAG",
    page_icon="📚",
    layout="wide"
)

# Initialize session state
if 'processed_files' not in st.session_state:
    st.session_state.processed_files = None
if 'knowledge_base' not in st.session_state:
    st.session_state.knowledge_base = []
if 'rag_processor' not in st.session_state:
    st.session_state.rag_processor = None

# Title and description
st.title("📚 Systematic Review Extractor with RAG")
st.markdown("""
This application helps you extract and analyze information from PDF files using AI and RAG (Retrieval-Augmented Generation).
Upload your PDFs to build a knowledge base and ask questions with enhanced context.
""")

# Get API key from secrets
api_key = st.secrets.get("OPENROUTER_API_KEY", "")
if not api_key:
    st.error("OpenRouter API key not found in secrets. Please add it to .streamlit/secrets.toml")
    st.stop()

class RAGProcessor:
    def __init__(self, api_key: str):
        self.api_key = api_key
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

# Initialize RAG processor
if st.session_state.rag_processor is None:
    st.session_state.rag_processor = RAGProcessor(api_key)

def process_pdf_directory(directory_path):
    """Process all PDF files in the given directory."""
    results = []
    for file_path in Path(directory_path).glob("*.pdf"):
        try:
            with open(file_path, 'rb') as file:
                # Here you would add your PDF processing logic
                # For now, we'll just return the filename
                results.append({
                    'filename': file_path.name,
                    'text': f"Sample text from {file_path.name}"
                })
        except Exception as e:
            st.error(f"Error processing {file_path.name}: {str(e)}")
    return results

def call_openrouter_api(prompt, api_key):
    """Call the OpenRouter API with the given prompt."""
    url = "https://openrouter.ai/api/v1/chat/completions"
    headers = {
        "Authorization": f"Bearer {api_key}",
        "HTTP-Referer": "https://wwwsystematicreviewextractor.streamlit.app",
        "Content-Type": "application/json"
    }
    payload = {
        "model": "meta-llama/llama-4-scout:free",
        "messages": [
            {"role": "system", "content": "You are a helpful AI assistant that extracts information from research papers."},
            {"role": "user", "content": prompt}
        ],
        "max_tokens": 4000
    }
    
    try:
        response = requests.post(url, headers=headers, json=payload)
        response.raise_for_status()
        return response.json()['choices'][0]['message']['content']
    except Exception as e:
        st.error(f"Error calling OpenRouter API: {str(e)}")
        return None

def process_text(text):
    """Process text using OpenRouter API."""
    prompt = f"""
    Extract the following information from this research paper text:
    - Title
    - Authors
    - Publication Year
    - Journal/Conference
    - Abstract
    - Key Findings
    - Methodology
    - Limitations
    
    Text: {text[:4000]}  # Limit text length to avoid token limits
    
    Respond with a JSON object containing these fields.
    """
    
    response = call_openrouter_api(prompt, api_key)
    if response:
        try:
            # Clean the response to ensure valid JSON
            response = response.strip()
            if response.startswith('```json'):
                response = response[7:]
            if response.endswith('```'):
                response = response[:-3]
            response = response.strip()
            
            return json.loads(response)
        except json.JSONDecodeError as e:
            st.error(f"Failed to parse AI response: {str(e)}")
            return {"error": str(e), "raw_response": response}
    return {"error": "No response from API"}

def analyze_findings(papers):
    """Analyze findings across multiple papers."""
    prompt = f"""
    Analyze these research papers and provide insights:
    {json.dumps(papers, indent=2)}
    
    Provide a structured analysis including:
    - Common themes
    - Research gaps
    - Key findings
    - Methodological patterns
    
    Respond with a JSON object containing these fields.
    """
    
    response = call_openrouter_api(prompt, api_key)
    if response:
        try:
            return json.loads(response)
        except json.JSONDecodeError as e:
            st.error(f"Failed to parse analysis response: {str(e)}")
            return {"error": str(e), "raw_response": response}
    return {"error": "No response from API"}

# Sidebar for RAG functionality
st.sidebar.title("RAG Knowledge Base")

# Knowledge base management
st.sidebar.subheader("Knowledge Base Management")

# Upload files for knowledge base
kb_files = st.sidebar.file_uploader(
    "Upload PDFs for Knowledge Base",
    type=ALLOWED_EXTENSIONS,
    accept_multiple_files=True,
    help="These files will be added to the RAG knowledge base"
)

if kb_files and st.sidebar.button("Add to Knowledge Base"):
    # Validate file sizes
    valid_kb_files = []
    for file in kb_files:
        if file.size > MAX_FILE_SIZE:
            st.sidebar.error(f"File {file.name} is too large. Maximum size is {MAX_FILE_SIZE/1024/1024}MB")
        else:
            valid_kb_files.append(file)
    
    if valid_kb_files:
        st.session_state.rag_processor.add_to_knowledge_base(valid_kb_files)

# Display knowledge base stats
stats = st.session_state.rag_processor.get_knowledge_base_stats()
st.sidebar.subheader("Knowledge Base Stats")
st.sidebar.write(f"Total Documents: {stats['total_documents']}")
st.sidebar.write(f"Total Chunks: {stats['total_chunks']}")
st.sidebar.write(f"Unique Files: {stats['unique_files']}")

# RAG Query Interface
st.sidebar.subheader("RAG Query")
rag_query = st.sidebar.text_area("Ask a question about your research papers:")
if st.sidebar.button("Ask RAG") and rag_query:
    if stats['total_chunks'] == 0:
        st.sidebar.error("Knowledge base is empty. Please add some PDFs first.")
    else:
        with st.spinner("Searching knowledge base..."):
            # Get relevant documents
            relevant_docs = st.session_state.rag_processor.search_knowledge_base(rag_query)
            
            if relevant_docs:
                st.sidebar.success(f"Found {len(relevant_docs)} relevant documents")
                
                # Display relevant documents
                with st.sidebar.expander("Relevant Documents"):
                    for i, doc in enumerate(relevant_docs):
                        st.write(f"**Document {i+1}:** {doc['metadata']['filename']}")
                        st.write(f"**Relevance:** {doc.get('similarity', 'N/A')}")
                        st.write(f"**Content:** {doc['text'][:200]}...")
                        st.write("---")
                
                # Generate RAG response
                with st.spinner("Generating RAG response..."):
                    enhanced_prompt = st.session_state.rag_processor.rag_query(rag_query)
                    rag_response = call_openrouter_api(enhanced_prompt, api_key)
                    
                    if rag_response:
                        st.sidebar.subheader("RAG Response")
                        st.sidebar.write(rag_response)
                    else:
                        st.sidebar.error("Failed to generate RAG response")
            else:
                st.sidebar.warning("No relevant documents found for your query")

# Main content area
st.subheader("PDF Processing")

# File uploader with size limit
uploaded_files = st.file_uploader(
    "Upload PDF files for processing",
    type=ALLOWED_EXTENSIONS,
    accept_multiple_files=True,
    help=f"Maximum file size: {MAX_FILE_SIZE/1024/1024}MB"
)

# Process uploaded files
if uploaded_files:
    # Validate file sizes
    valid_files = []
    for file in uploaded_files:
        if file.size > MAX_FILE_SIZE:
            st.error(f"File {file.name} is too large. Maximum size is {MAX_FILE_SIZE/1024/1024}MB")
        else:
            valid_files.append(file)
    
    if valid_files:
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Save uploaded files
            for file in valid_files:
                with open(temp_path / file.name, "wb") as f:
                    f.write(file.getvalue())
            
            if st.button("Process PDFs"):
                with st.spinner("Processing PDFs..."):
                    try:
                        # Process PDFs
                        results = process_pdf_directory(str(temp_path))
                        
                        if results:
                            # Convert results to DataFrame
                            df = pd.DataFrame(results)
                            
                            # Process with AI
                            st.subheader("AI Analysis")
                            
                            # Create progress bar
                            progress_bar = st.progress(0)
                            status_text = st.empty()
                            
                            # Process each PDF with AI
                            ai_results = []
                            total_files = len(df)
                            
                            for idx, (_, row) in enumerate(df.iterrows()):
                                status_text.text(f"Processing file {idx + 1} of {total_files}")
                                progress_bar.progress((idx + 1) / total_files)
                                
                                text = row['text'] if 'text' in row else ""
                                if text:
                                    st.write(f"Processing: {row.get('filename', 'unknown')}")
                                    ai_data = process_text(text)
                                    if "error" in ai_data:
                                        st.error(f"Error processing file {row.get('filename', 'unknown')}:")
                                        st.error(f"Error details: {ai_data['error']}")
                                        if "raw_response" in ai_data:
                                            with st.expander("View raw response"):
                                                st.code(ai_data['raw_response'])
                                    ai_results.append(ai_data)
                            
                            # Clear progress indicators
                            progress_bar.empty()
                            status_text.empty()
                            
                            # Analyze findings across all papers
                            if ai_results:
                                with st.spinner("Analyzing findings..."):
                                    analysis = analyze_findings(ai_results)
                                    
                                    # Display analysis
                                    st.markdown("### Research Analysis")
                                    if isinstance(analysis.get('analysis', {}), dict):
                                        for key, value in analysis['analysis'].items():
                                            st.markdown(f"#### {key.replace('_', ' ').title()}")
                                            if isinstance(value, list):
                                                for item in value:
                                                    st.markdown(f"- {item}")
                                            else:
                                                st.write(value)
                                    else:
                                        st.write(analysis.get('analysis', 'No analysis available'))
                                    
                                    # Create AI-enhanced DataFrame
                                    ai_df = pd.DataFrame(ai_results)
                                    
                                    # Save both DataFrames to a temporary file
                                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                                    excel_filename = f"pdf_extracts_ai_{timestamp}.xlsx"
                                    
                                    with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as tmp:
                                        with pd.ExcelWriter(tmp.name) as writer:
                                            df.to_excel(writer, sheet_name='Raw Data', index=False)
                                            ai_df.to_excel(writer, sheet_name='AI Analysis', index=False)
                                        
                                        # Store the file path in session state
                                        st.session_state.processed_files = {
                                            'path': tmp.name,
                                            'filename': excel_filename
                                        }
                                    
                                    st.success("Files processed successfully!")
                                    
                                    # Display preview of AI analysis
                                    st.subheader("Preview of AI Analysis")
                                    st.dataframe(ai_df)
                                    
                                    # Add download button
                                    with open(st.session_state.processed_files['path'], 'rb') as f:
                                        st.download_button(
                                            label="Download Excel File",
                                            data=f,
                                            file_name=st.session_state.processed_files['filename'],
                                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                        )
                            else:
                                st.error("No text content found in PDFs for AI analysis")
                        else:
                            st.error("No results found from PDF processing")
                    except Exception as e:
                        st.error(f"Error processing PDFs: {str(e)}")
                        st.exception(e)

# Add information about the app
st.sidebar.markdown("""
### About
This app uses AI and RAG to:
- Extract structured information from PDFs
- Build a searchable knowledge base
- Answer questions with context
- Analyze research findings
- Identify patterns and insights

### Requirements
- PDF files to analyze (max 10MB each)
- OpenRouter API key
""")

# Clean up temporary files on session end
if st.session_state.processed_files and os.path.exists(st.session_state.processed_files['path']):
    try:
        os.unlink(st.session_state.processed_files['path'])
    except:
        pass 