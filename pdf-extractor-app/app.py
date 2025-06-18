import os
import sys
import streamlit as st

# Add the project root to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from src.agent.pdf_extractor_streamlit import *
except ImportError as e:
    st.error(f"Error importing the application: {str(e)}")
    st.error("Please make sure all required files are present in the correct locations.")
    st.stop()

# This file serves as the entry point for Streamlit Cloud
# It simply imports and runs the main app 