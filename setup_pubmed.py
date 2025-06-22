#!/usr/bin/env python3
"""
Setup script for PubMed integration
This script helps configure the email address required for PubMed access.
"""

import os
import getpass

def setup_pubmed_email():
    """Set up email for PubMed access."""
    print("🔧 PubMed Setup")
    print("=" * 40)
    print("NCBI requires an email address for PubMed access.")
    print("This helps them contact you if there are issues with your usage.")
    print()
    
    # Check if email is already set
    current_email = os.environ.get("PUBMED_EMAIL")
    if current_email:
        print(f"Current email: {current_email}")
        change = input("Do you want to change it? (y/n): ").lower()
        if change != 'y':
            return current_email
    
    # Get email from user
    while True:
        email = input("Enter your email address: ").strip()
        if '@' in email and '.' in email:
            break
        print("❌ Please enter a valid email address.")
    
    # Save to environment variable
    os.environ["PUBMED_EMAIL"] = email
    
    # Update the main script
    update_main_script(email)
    
    print(f"✅ Email configured: {email}")
    print("You can now use PubMed integration!")
    return email

def update_main_script(email):
    """Update the main script with the email address."""
    try:
        with open("run_rag_interactive.py", "r") as f:
            content = f.read()
        
        # Replace the placeholder email
        updated_content = content.replace(
            'Entrez.email = "your-email@example.com"',
            f'Entrez.email = "{email}"'
        )
        
        with open("run_rag_interactive.py", "w") as f:
            f.write(updated_content)
            
        print("✅ Main script updated with your email address.")
        
    except Exception as e:
        print(f"❌ Error updating script: {e}")
        print("Please manually update the email in run_rag_interactive.py")

if __name__ == "__main__":
    setup_pubmed_email() 