import streamlit as st
import os
from pathlib import Path
import pandas as pd
from datetime import datetime
from math_example import process_pdf_directory
import time

st.set_page_config(
    page_title="Systematic Review PDF Extractor",
    page_icon="📚",
    layout="wide"
)

def main():
    st.title("📚 Systematic Review PDF Extractor")
    st.markdown("""
    This application helps you extract information from PDF files for systematic reviews.
    Upload your PDF files and protocol, and the app will process them and generate an Excel file with the extracted information.
    """)

    # Create two columns for file uploads
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("PDF Files")
        pdf_files = st.file_uploader(
            "Upload your PDF files",
            type=['pdf'],
            accept_multiple_files=True
        )

    with col2:
        st.subheader("Protocol File")
        protocol_file = st.file_uploader(
            "Upload your protocol file",
            type=['pdf'],
            accept_multiple_files=False
        )

    # Output directory selection
    st.subheader("Output Settings")
    output_dir = st.text_input(
        "Output Directory",
        value=str(Path.home() / "Downloads"),
        help="Where to save the Excel file"
    )

    # Process button
    if st.button("Start Extraction", type="primary"):
        if not pdf_files:
            st.error("Please upload at least one PDF file")
            return
        
        if not protocol_file:
            st.error("Please upload a protocol file")
            return

        # Create a temporary directory for uploaded files
        temp_dir = Path("temp_uploads")
        temp_dir.mkdir(exist_ok=True)

        try:
            # Save uploaded files
            pdf_paths = []
            for pdf_file in pdf_files:
                pdf_path = temp_dir / pdf_file.name
                with open(pdf_path, "wb") as f:
                    f.write(pdf_file.getbuffer())
                pdf_paths.append(pdf_path)

            protocol_path = temp_dir / protocol_file.name
            with open(protocol_path, "wb") as f:
                f.write(protocol_file.getbuffer())

            # Create progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()

            # Process files
            status_text.text("Processing files...")
            output_file = process_pdf_directory(str(temp_dir), str(protocol_path))

            if output_file and os.path.exists(output_file):
                progress_bar.progress(100)
                status_text.text("Extraction completed successfully!")
                
                # Show success message and download button
                st.success(f"✅ Extraction completed! File saved to: {output_file}")
                
                # Create download button
                with open(output_file, "rb") as f:
                    st.download_button(
                        label="Download Excel File",
                        data=f,
                        file_name=os.path.basename(output_file),
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
            else:
                st.error("Extraction failed - no output file was created")

        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
        
        finally:
            # Clean up temporary files
            for file in temp_dir.glob("*"):
                file.unlink()
            temp_dir.rmdir()

if __name__ == "__main__":
    main() 