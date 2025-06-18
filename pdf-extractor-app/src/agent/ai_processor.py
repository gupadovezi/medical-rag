import os
from typing import List, Dict, Any
from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain.output_parsers import PydanticOutputParser
from pydantic import BaseModel, Field
import json

# Load environment variables
load_dotenv()

class ExtractedData(BaseModel):
    """Schema for extracted data from PDFs"""
    title: str = Field(description="Title of the paper")
    authors: List[str] = Field(description="List of authors")
    year: int = Field(description="Publication year")
    abstract: str = Field(description="Abstract of the paper")
    key_findings: List[str] = Field(description="Key findings or conclusions")
    methodology: str = Field(description="Research methodology used")
    keywords: List[str] = Field(description="Keywords or topics")

class AIProcessor:
    def __init__(self):
        """Initialize the AI processor with OpenAI model"""
        self.llm = ChatOpenAI(
            model="gpt-3.5-turbo",
            temperature=0,
            api_key=os.getenv("OPENAI_API_KEY")
        )
        self.parser = PydanticOutputParser(pydantic_object=ExtractedData)
        
        # Create the prompt template
        self.prompt = ChatPromptTemplate.from_messages([
            ("system", """You are an expert research paper analyzer. 
            Extract key information from the provided text following this format:
            {format_instructions}
            
            Focus on extracting accurate information and maintain objectivity."""),
            ("user", "{text}")
        ])

    def process_text(self, text: str) -> Dict[str, Any]:
        """Process text using AI to extract structured information"""
        try:
            # Format the prompt with instructions and text
            formatted_prompt = self.prompt.format_messages(
                format_instructions=self.parser.get_format_instructions(),
                text=text
            )
            
            # Get response from the model
            response = self.llm.invoke(formatted_prompt)
            
            # Parse the response
            parsed_data = self.parser.parse(response.content)
            
            return parsed_data.dict()
            
        except Exception as e:
            print(f"Error processing text with AI: {str(e)}")
            return {
                "error": str(e),
                "title": "",
                "authors": [],
                "year": None,
                "abstract": "",
                "key_findings": [],
                "methodology": "",
                "keywords": []
            }

    def analyze_findings(self, extracted_data: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze multiple papers to find patterns and insights"""
        try:
            # Create a summary prompt
            summary_prompt = ChatPromptTemplate.from_messages([
                ("system", """You are an expert research analyst. 
                Analyze the following research papers and provide:
                1. Common themes
                2. Conflicting findings
                3. Research gaps
                4. Future research directions
                
                Be objective and data-driven in your analysis."""),
                ("user", "{papers}")
            ])
            
            # Format papers data
            papers_text = json.dumps(extracted_data, indent=2)
            
            # Get analysis from the model
            formatted_prompt = summary_prompt.format_messages(papers=papers_text)
            response = self.llm.invoke(formatted_prompt)
            
            return {
                "analysis": response.content,
                "papers_analyzed": len(extracted_data)
            }
            
        except Exception as e:
            print(f"Error analyzing findings: {str(e)}")
            return {
                "error": str(e),
                "analysis": "",
                "papers_analyzed": 0
            } 