import streamlit as st
import os
from pathlib import Path
import pandas as pd
from datetime import datetime
from math_example import process_pdf_directory
from ai_processor import AIProcessor

st.set_page_config(
    page_title="AI-Powered PDF Extractor",
    page_icon="📚",
    layout="wide"
)

# Initialize session state
if 'ai_processor' not in st.session_state:
    st.session_state.ai_processor = None

def main():
    st.title("📚 AI-Powered PDF Extractor")
    st.markdown("""
    This application helps you extract and analyze information from PDF files using AI.
    Upload your PDFs and let the AI process them to extract key information and insights.
    """)

    # Add API key configuration
    st.sidebar.title("Configuration")
    api_key = st.sidebar.text_input(
        "OpenAI API Key",
        value="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",
        type="password"
    )

    if api_key:
        os.environ["OPENAI_API_KEY"] = api_key
        if st.session_state.ai_processor is None:
            st.session_state.ai_processor = AIProcessor()
        st.sidebar.success("API key configured!")
    else:
        st.sidebar.warning("Please enter your OpenAI API key to enable AI features")

    # File uploader
    uploaded_files = st.file_uploader("Upload PDF files", type=['pdf'], accept_multiple_files=True)

    # Output directory selection
    output_dir = st.text_input(
        "Output Directory",
        value=str(Path.home() / "Downloads"),
        help="Directory where the Excel file will be saved"
    )

    if uploaded_files:
        # Create a temporary directory for uploaded files
        temp_dir = Path("temp_uploads")
        temp_dir.mkdir(exist_ok=True)
        
        # Save uploaded files
        for uploaded_file in uploaded_files:
            with open(temp_dir / uploaded_file.name, "wb") as f:
                f.write(uploaded_file.getvalue())
        
        if st.button("Process PDFs"):
            if not st.session_state.ai_processor:
                st.error("Please configure your API key in the sidebar first")
            else:
                with st.spinner("Processing PDFs..."):
                    try:
                        # Process PDFs without protocol
                        results = process_pdf_directory(str(temp_dir))
                        
                        if results:
                            # Convert results to DataFrame
                            df = pd.DataFrame(results)
                            
                            # Process with AI
                            st.subheader("AI Analysis")
                            with st.spinner("Analyzing with AI..."):
                                # Process each PDF with AI
                                ai_results = []
                                for _, row in df.iterrows():
                                    text = row['text'] if 'text' in row else ""
                                    if text:
                                        ai_data = st.session_state.ai_processor.process_text(text)
                                        ai_results.append(ai_data)
                                
                                # Analyze findings across all papers
                                if ai_results:
                                    analysis = st.session_state.ai_processor.analyze_findings(ai_results)
                                    
                                    # Display analysis
                                    st.markdown("### Research Analysis")
                                    st.write(analysis['analysis'])
                                    
                                    # Create AI-enhanced DataFrame
                                    ai_df = pd.DataFrame(ai_results)
                                    
                                    # Save both DataFrames
                                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                                    excel_path = Path(output_dir) / f"pdf_extracts_ai_{timestamp}.xlsx"
                                    
                                    with pd.ExcelWriter(excel_path) as writer:
                                        df.to_excel(writer, sheet_name='Raw Data', index=False)
                                        ai_df.to_excel(writer, sheet_name='AI Analysis', index=False)
                                    
                                    st.success(f"Files processed successfully! Saved to: {excel_path}")
                                    
                                    # Display preview of AI analysis
                                    st.subheader("Preview of AI Analysis")
                                    st.dataframe(ai_df)
                                else:
                                    st.error("No text content found in PDFs for AI analysis")
                        else:
                            st.error("No results found from PDF processing")
                    except Exception as e:
                        st.error(f"Error processing PDFs: {str(e)}")
                
                # Clean up temporary files
                for file in temp_dir.glob("*.pdf"):
                    file.unlink()
                temp_dir.rmdir()

    # Add information about the app
    st.sidebar.markdown("""
    ### About
    This app uses AI to:
    - Extract structured information from PDFs
    - Analyze research findings
    - Identify patterns and insights
    - Generate comprehensive reports

    ### Requirements
    - OpenAI API key
    - PDF files to analyze
    """)

if __name__ == "__main__":
    main() 