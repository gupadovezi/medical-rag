import streamlit as st
import pandas as pd
from pathlib import Path
import os
from datetime import datetime
import tempfile
import shutil
from math_example import process_pdf_directory
from ai_processor import AIProcessor

# Constants
MAX_FILE_SIZE = 10 * 1024 * 1024  # 10MB
ALLOWED_EXTENSIONS = {'pdf'}

# Set page config
st.set_page_config(
    page_title="AI-Powered PDF Extractor",
    page_icon="📚",
    layout="wide"
)

# Initialize session state
if 'ai_processor' not in st.session_state:
    st.session_state.ai_processor = None
if 'processed_files' not in st.session_state:
    st.session_state.processed_files = None

# Title and description
st.title("📚 AI-Powered PDF Extractor")
st.markdown("""
This application helps you extract and analyze information from PDF files using AI.
Upload your PDFs and let the AI process them to extract key information and insights.
""")

# Get API key from secrets
api_key = st.secrets.get("OPENROUTER_API_KEY", "")
if not api_key:
    st.error("OpenRouter API key not found in secrets. Please add it to .streamlit/secrets.toml")
    st.stop()

# Initialize AI processor
try:
    os.environ["OPENROUTER_API_KEY"] = api_key
    st.session_state.ai_processor = AIProcessor()
except Exception as e:
    st.error(f"Error initializing AI processor: {str(e)}")
    st.stop()

# File uploader with size limit
uploaded_files = st.file_uploader(
    "Upload PDF files",
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
                                    ai_data = st.session_state.ai_processor.process_text(text)
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
                                    analysis = st.session_state.ai_processor.analyze_findings(ai_results)
                                    
                                    # Display analysis
                                    st.markdown("### Research Analysis")
                                    if isinstance(analysis['analysis'], dict):
                                        for key, value in analysis['analysis'].items():
                                            st.markdown(f"#### {key.replace('_', ' ').title()}")
                                            if isinstance(value, list):
                                                for item in value:
                                                    st.markdown(f"- {item}")
                                            else:
                                                st.write(value)
                                    else:
                                        st.write(analysis['analysis'])
                                    
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
This app uses AI to:
- Extract structured information from PDFs
- Analyze research findings
- Identify patterns and insights
- Generate comprehensive reports

### Requirements
- PDF files to analyze (max 10MB each)
""")

# Clean up temporary files on session end
if st.session_state.processed_files and os.path.exists(st.session_state.processed_files['path']):
    try:
        os.unlink(st.session_state.processed_files['path'])
    except:
        pass 