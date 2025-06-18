import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from pathlib import Path
import threading
from math_example import process_pdf_directory
from datetime import datetime

class PDFExtractorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Systematic Review PDF Extractor")
        self.root.geometry("800x600")
        
        # Create main frame
        self.main_frame = ttk.Frame(root, padding="10")
        self.main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # PDF Directory Selection
        ttk.Label(self.main_frame, text="PDF Directory:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.pdf_dir_var = tk.StringVar()
        self.pdf_dir_entry = ttk.Entry(self.main_frame, textvariable=self.pdf_dir_var, width=50)
        self.pdf_dir_entry.grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(self.main_frame, text="Browse", command=self.browse_pdf_dir).grid(row=0, column=2, padx=5)
        
        # Protocol File Selection
        ttk.Label(self.main_frame, text="Protocol File:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.protocol_var = tk.StringVar()
        self.protocol_entry = ttk.Entry(self.main_frame, textvariable=self.protocol_var, width=50)
        self.protocol_entry.grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(self.main_frame, text="Browse", command=self.browse_protocol).grid(row=1, column=2, padx=5)
        
        # Output Directory Selection
        ttk.Label(self.main_frame, text="Output Directory:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.output_dir_var = tk.StringVar()
        self.output_dir_entry = ttk.Entry(self.main_frame, textvariable=self.output_dir_var, width=50)
        self.output_dir_entry.grid(row=2, column=1, padx=5, pady=5)
        ttk.Button(self.main_frame, text="Browse", command=self.browse_output_dir).grid(row=2, column=2, padx=5)
        
        # Progress Frame
        self.progress_frame = ttk.LabelFrame(self.main_frame, text="Progress", padding="5")
        self.progress_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)
        
        # Progress Bar
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(self.progress_frame, variable=self.progress_var, maximum=100)
        self.progress_bar.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=5, pady=5)
        
        # Status Label
        self.status_var = tk.StringVar(value="Ready")
        self.status_label = ttk.Label(self.progress_frame, textvariable=self.status_var)
        self.status_label.grid(row=1, column=0, sticky=tk.W, padx=5)
        
        # Output File Label
        self.output_file_var = tk.StringVar(value="Output will be saved to: Not selected")
        self.output_file_label = ttk.Label(self.progress_frame, textvariable=self.output_file_var)
        self.output_file_label.grid(row=2, column=0, sticky=tk.W, padx=5)
        
        # Log Frame
        self.log_frame = ttk.LabelFrame(self.main_frame, text="Log", padding="5")
        self.log_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
        
        # Log Text
        self.log_text = tk.Text(self.log_frame, height=15, width=80)
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Scrollbar for Log
        scrollbar = ttk.Scrollbar(self.log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.log_text['yscrollcommand'] = scrollbar.set
        
        # Start Button
        self.start_button = ttk.Button(self.main_frame, text="Start Extraction", command=self.start_extraction)
        self.start_button.grid(row=5, column=0, columnspan=3, pady=10)
        
        # Configure grid weights
        self.main_frame.columnconfigure(1, weight=1)
        self.main_frame.rowconfigure(4, weight=1)
        self.log_frame.columnconfigure(0, weight=1)
        self.log_frame.rowconfigure(0, weight=1)
        
        # Set default output directory to user's Downloads folder
        default_output = str(Path.home() / "Downloads")
        self.output_dir_var.set(default_output)
        self.update_output_file_label()
        
    def browse_pdf_dir(self):
        directory = filedialog.askdirectory()
        if directory:
            self.pdf_dir_var.set(directory)
            self.log("Selected PDF directory: " + directory)
    
    def browse_protocol(self):
        filename = filedialog.askopenfilename(
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")]
        )
        if filename:
            self.protocol_var.set(filename)
            self.log("Selected protocol file: " + filename)
    
    def browse_output_dir(self):
        directory = filedialog.askdirectory()
        if directory:
            self.output_dir_var.set(directory)
            self.log("Selected output directory: " + directory)
            self.update_output_file_label()
    
    def update_output_file_label(self):
        output_dir = self.output_dir_var.get()
        if output_dir:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = os.path.join(output_dir, f"systematic_review_extracts_{timestamp}.xlsx")
            self.output_file_var.set(f"Output will be saved to: {output_file}")
    
    def log(self, message):
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
    
    def update_progress(self, value, message):
        self.progress_var.set(value)
        self.status_var.set(message)
        self.log(message)
    
    def start_extraction(self):
        pdf_dir = self.pdf_dir_var.get()
        protocol_path = self.protocol_var.get()
        output_dir = self.output_dir_var.get()
        
        if not pdf_dir or not protocol_path:
            messagebox.showerror("Error", "Please select both PDF directory and protocol file")
            return
        
        if not os.path.exists(pdf_dir):
            messagebox.showerror("Error", "PDF directory does not exist")
            return
        
        if not os.path.exists(protocol_path):
            messagebox.showerror("Error", "Protocol file does not exist")
            return
        
        # Create output directory if it doesn't exist
        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Error", f"Could not create output directory: {str(e)}")
            return
        
        # Disable start button during processing
        self.start_button.state(['disabled'])
        
        # Start extraction in a separate thread
        thread = threading.Thread(target=self.run_extraction, args=(pdf_dir, protocol_path, output_dir))
        thread.daemon = True
        thread.start()
    
    def run_extraction(self, pdf_dir, protocol_path, output_dir):
        try:
            self.update_progress(0, "Starting extraction...")
            
            # Get list of PDF files
            pdf_files = list(Path(pdf_dir).glob("*.pdf"))
            total_files = len(pdf_files)
            
            if total_files == 0:
                self.update_progress(0, "No PDF files found in the selected directory")
                return
            
            self.log(f"Found {total_files} PDF files")
            
            # Process files
            for i, pdf_file in enumerate(pdf_files, 1):
                progress = (i / total_files) * 100
                self.update_progress(progress, f"Processing {pdf_file.name} ({i}/{total_files})")
            
            # Run the extraction
            output_file = process_pdf_directory(pdf_dir, protocol_path)
            
            if output_file and os.path.exists(output_file):
                self.update_progress(100, "Extraction completed successfully!")
                self.log(f"Output saved to: {output_file}")
                messagebox.showinfo("Success", f"PDF extraction completed successfully!\nOutput saved to: {output_file}")
            else:
                self.update_progress(0, "Extraction failed - no output file was created")
                self.log("Error: No output file was created")
                messagebox.showerror("Error", "Extraction failed - no output file was created")
            
        except Exception as e:
            self.log(f"Error: {str(e)}")
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
        
        finally:
            # Re-enable start button
            self.start_button.state(['!disabled'])

def main():
    root = tk.Tk()
    app = PDFExtractorGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main() 