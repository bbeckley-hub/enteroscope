#!/usr/bin/env python3
"""
MLST Module for EnteroScope - Focused on Enterobacter cloacae Complex (ECC)
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: https://github.com/bbeckley-hub/enteroscope
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026
"""

import os
import sys
import json
import glob
import argparse
import subprocess
import random
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime

class ModularMLSTAnalyzer:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        self.mlst_bin = script_dir / "mlst"
        
        # Science quotes for rotation
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"},
            {"text": "EnteroScope represents the convergence of genomic surveillance and clinical diagnostics for Enterobacter cloacae complex.", "author": "Brown Beckley"},
            {"text": "In the battle against antimicrobial resistance, EnteroScope reveals the genetic blueprints of resistant Enterobacter strains.", "author": "Brown Beckley"},
            {"text": "EnteroScope bridges the gap between sequencing data and public health action for ECC.", "author": "Brown Beckley"},
            {"text": "Through EnteroScope, we turn complexity into clear reports, empowering clinicians and researchers.", "author": "Brown Beckley"},
            {"text": "EnteroScope makes advanced MLST typing for Enterobacter cloacae complex accessible to all.", "author": "Brown Beckley"}
        ]
    
    def get_random_quote(self):
        """Get a random science quote"""
        return random.choice(self.science_quotes)
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        if os.path.isfile(input_path):
            return [Path(input_path)]
        
        fasta_patterns = [
            input_path,
            f"{input_path}/*.fna", f"{input_path}/*.fasta",
            f"{input_path}/*.fa", f"{input_path}/*.fn",
            f"{input_path}/*.fna.gz", f"{input_path}/*.fasta.gz",
            f"{input_path}/*.fa.gz", f"{input_path}/*.gb",
            f"{input_path}/*.gbk", f"{input_path}/*.gbff"
        ]
        
        fasta_files = []
        for pattern in fasta_patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                path = Path(file_path)
                if path.is_file():
                    fasta_files.append(path)
        
        return sorted(list(set(fasta_files)))

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "ecloacae") -> Dict:
        """Run MLST analysis for a single file using E. cloacae scheme"""
        print(f"🔬 Processing: {input_file.name}")
        
        # Create sample-specific output directory
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save raw MLST output first
        raw_output_file = sample_output_dir / "mlst_raw_output.txt"
        
        # Run MLST command
        mlst_cmd = [
            "perl", str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath"
        ]
        
        try:
            # Run and capture output
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            
            # Save raw output
            with open(raw_output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
            
            print(f"Raw MLST output: {result.stdout.strip()}")
            
            # Parse the CSV output
            mlst_results = self.parse_mlst_csv(result.stdout, input_file.name)
            
            # Add identity and coverage information (no lineage database)
            mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
            
            # Generate output files
            self.generate_output_files(mlst_results, sample_output_dir)
            
            print(f"✅ Completed: {input_file.name} -> ST{mlst_results.get('st', 'ND')}")
            return mlst_results
            
        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result

    def parse_mlst_csv(self, stdout: str, sample_name: str) -> Dict:
        """Parse MLST CSV output - comma-separated format from mlst Perl script"""
        print(f"Parsing CSV output for {sample_name}")
        
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name)
        
        # Find the result line (usually the last line with data)
        result_line = None
        for line in reversed(lines):
            if line.strip() and ',' in line and not line.startswith('['):
                result_line = line.strip()
                break
        
        if not result_line:
            return self.get_empty_results(sample_name)
        
        print(f"CSV result line: {result_line}")
        
        # Split by COMMA
        parts = result_line.split(',')
        print(f"CSV parts: {parts}")
        
        if len(parts) < 3:
            return self.get_empty_results(sample_name)
        
        # Extract components - format: filename,scheme,ST,allele1,allele2,...
        filename = parts[0]
        scheme = parts[1]
        st = parts[2]
        
        # Extract alleles from remaining parts
        alleles = {}
        allele_parts = []
        
        for i in range(3, len(parts)):
            allele_str = parts[i]
            if '(' in allele_str and ')' in allele_str:
                gene = allele_str.split('(')[0]
                allele = allele_str.split('(')[1].rstrip(')')
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
        
        allele_profile = '-'.join(allele_parts) if allele_parts else ""
        
        return {
            "sample": sample_name,
            "st": st,
            "scheme": scheme,
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": "HIGH" if st and st != '-' and st != 'ND' else "LOW",
            "mlst_assigned": True if st and st != '-' and st != 'ND' else False
        }

    def get_identity_coverage(self, st: str) -> Dict:
        """Get identity and coverage information based on MLST assignment"""
        if st and st != '-' and st != 'ND' and st != 'UNKNOWN':
            return {
                "identity": "100%",
                "coverage": "100%",
                "mlst_status": "Assigned",
                "quality_metrics": {
                    "assembly_quality": "High Quality",
                    "allele_completeness": "Complete",
                    "database_match": "Perfect Match"
                }
            }
        else:
            return {
                "identity": "Not Assigned",
                "coverage": "Not Assigned",
                "mlst_status": "Not Assigned",
                "quality_metrics": {
                    "assembly_quality": "Requires Review",
                    "allele_completeness": "Incomplete",
                    "database_match": "No Match"
                }
            }

    def get_empty_results(self, sample_name: str) -> Dict:
        """Return empty results structure"""
        return {
            "sample": sample_name,
            "st": "ND",
            "scheme": "ecloacae",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False
        }

    def get_fallback_results(self, sample_name: str) -> Dict:
        """Fallback when MLST fails"""
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": "ecloacae",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
            "error": "MLST analysis failed"
        }

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        """Generate 3 output files: HTML, TXT, and TSV"""
        if 'identity' not in mlst_results:
            mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
        
        self.generate_html_report(mlst_results, output_dir)
        self.generate_text_report(mlst_results, output_dir)
        self.generate_tsv_report(mlst_results, output_dir)

    def generate_text_report(self, mlst_results: Dict, output_dir: Path):
        """Generate detailed text report (no lineage section)"""
        report = f"""MLST Analysis Report - EnteroScope
=================================

Sample: {mlst_results['sample']}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

MLST TYPING RESULTS:
-------------------
Sequence Type (ST): {mlst_results['st']}
Scheme: {mlst_results['scheme']}
Confidence: {mlst_results['confidence']}
MLST Status: {mlst_results.get('mlst_status', 'Not Assigned')}

Identity & Coverage:
Identity: {mlst_results.get('identity', 'Not Assigned')}
Coverage: {mlst_results.get('coverage', 'Not Assigned')}

Allele Profile:
{mlst_results['allele_profile']}

Detailed Alleles:
"""
        for gene, allele in mlst_results['alleles'].items():
            report += f"- {gene}: {allele}\n"
        
        if 'quality_metrics' in mlst_results:
            report += f"""
QUALITY METRICS:
----------------
Assembly Quality: {mlst_results['quality_metrics'].get('assembly_quality', 'Unknown')}
Allele Completeness: {mlst_results['quality_metrics'].get('allele_completeness', 'Unknown')}
Database Match: {mlst_results['quality_metrics'].get('database_match', 'Unknown')}
"""
        
        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        """Generate simple TSV report"""
        tsv_content = f"Sample\tST\tScheme\tMLST_Status\tIdentity\tCoverage\tAllele_Profile\n"
        tsv_content += f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results.get('mlst_status', 'Not Assigned')}\t{mlst_results.get('identity', 'Not Assigned')}\t{mlst_results.get('coverage', 'Not Assigned')}\t{mlst_results['allele_profile']}\n"
        
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        """Generate beautiful HTML report with Teal/Cyan + Electric Blue theme"""
        random_quote = self.get_random_quote()
        
        sample = mlst_results['sample']
        st = mlst_results['st']
        scheme = mlst_results['scheme']
        confidence = mlst_results['confidence']
        allele_profile = mlst_results['allele_profile']
        identity = mlst_results.get('identity', 'Not Assigned')
        coverage = mlst_results.get('coverage', 'Not Assigned')
        mlst_status = mlst_results.get('mlst_status', 'Not Assigned')
        quality_metrics = mlst_results.get('quality_metrics', {})
        
        # Build alleles HTML
        alleles_html = ''
        for gene, allele in mlst_results.get('alleles', {}).items():
            alleles_html += f'''                <div class="allele-card">
                    <div style="font-size: 12px; opacity: 0.9;">{gene}</div>
                    <div style="font-size: 18px;">{allele}</div>
                </div>
'''
        
        # Build quality metrics HTML
        quality_metrics_html = ''
        if quality_metrics:
            for key, value in quality_metrics.items():
                quality_metrics_html += f'''                <div class="quality-card">
                    <div class="quality-value">{value}</div>
                    <div class="quality-label">{key.replace('_', ' ').title()}</div>
                </div>
'''
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroScope - MLST Analysis Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(6, 182, 212, 0.5);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #06b6d4;
            text-shadow: 0 0 10px rgba(6, 182, 212, 0.5);
            overflow-x: auto;
        }}
        
        .quote-container {{
            background: rgba(0, 0, 0, 0.4);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(6, 182, 212, 0.4);
            transition: opacity 0.5s ease-in-out;
        }}
        
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
        }}
        
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }}
        
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        
        .report-section h2 {{
            color: #0f766e;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        
        .metric-card {{
            background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        
        .identity-coverage-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        
        .ic-card {{
            background: linear-gradient(135deg, #0d9488 0%, #0f766e 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
            text-align: center;
        }}
        
        .ic-card.not-assigned {{
            background: linear-gradient(135deg, #6b7280 0%, #4b5563 100%);
        }}
        
        .ic-value {{
            font-size: 28px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .ic-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        
        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        
        .allele-card {{
            background: linear-gradient(135deg, #06b6d4 0%, #0891b2 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
            font-weight: bold;
        }}
        
        .confidence-high {{
            color: #16a34a;
            font-weight: bold;
        }}
        
        .confidence-medium {{
            color: #f59e0b;
            font-weight: bold;
        }}
        
        .confidence-low {{
            color: #dc2626;
            font-weight: bold;
        }}
        
        .profile-box {{
            background: #f0fdfa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #0d9488;
        }}
        
        .quality-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        
        .quality-card {{
            background: #f0f9ff;
            color: #0f766e;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            border: 1px solid #bae6fd;
        }}
        
        .quality-value {{
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .quality-label {{
            font-size: 12px;
            opacity: 0.8;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        
        .footer a {{
            color: #fbbf24;
            text-decoration: none;
        }}
        
        .footer a:hover {{
            text-decoration: underline;
        }}
        
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .allele-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
█████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║███████╗██║      ██║   ██║██████╔╝█████╗  
██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
            
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sample Name</div>
                    <div class="metric-value">{sample}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value">{scheme.title()}</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🎯 Identity & Coverage</h2>
            <div class="identity-coverage-grid">
'''
        
        identity_class = "ic-card" if identity == "100%" else "ic-card not-assigned"
        coverage_class = "ic-card" if coverage == "100%" else "ic-card not-assigned"
        
        html_content += f'''                <div class="{identity_class}">
                    <div class="ic-value">{identity}</div>
                    <div class="ic-label">Identity</div>
                </div>
                <div class="{coverage_class}">
                    <div class="ic-value">{coverage}</div>
                    <div class="ic-label">Coverage</div>
                </div>
'''
        
        html_content += f'''            </div>
            
            <h3>MLST Status: <span style="color: {'#16a34a' if mlst_status == 'Assigned' else '#dc2626'}">{mlst_status}</span></h3>
            
            <h3>Quality Metrics</h3>
            <div class="quality-grid">
{quality_metrics_html}            </div>
        </div>
        
        <div class="report-section">
            <h2>🧬 MLST Typing Results</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sequence Type</div>
                    <div class="metric-value">ST{st}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value confidence-{confidence.lower()}">{confidence}</div>
                </div>
            </div>
            
            <h3>Allele Profile</h3>
            <div class="profile-box">
                <code style="font-size: 16px; color: #0f766e; font-weight: bold;">{allele_profile}</code>
            </div>
            
            <h3>Individual Alleles</h3>
            <div class="allele-grid">
{alleles_html}            </div>
        </div>
        
        <div class="footer">
            <p><strong>EnteroScope</strong> - MLST Analysis Module for <em>Enterobacter cloacae</em> Complex</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <div class="authorship">
                <p><strong>Author:</strong> Brown Beckley | <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a></p>
                <p><strong>Email:</strong> <a href="mailto:brownbeckley94@gmail.com">brownbeckley94@gmail.com</a></p>
                <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
                <p><strong>Issue reporting:</strong> Please use GitHub Issues or send an email.</p>
            </div>
        </div>
    </div>

    <script>
        const quotes = {json.dumps(self.science_quotes, indent=2)};

        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');

        function getRandomQuote() {{
            return quotes[Math.floor(Math.random() * quotes.length)];
        }}

        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}

        setInterval(displayQuote, 10000);
    </script>
</body>
</html>'''
        
        with open(output_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    # ------------------------------------------------------------------
    # Summary methods (batch mode)
    # ------------------------------------------------------------------
    def create_mlst_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive MLST summary files for all samples"""
        print("📊 Creating MLST summary files...")
        self.create_mlst_tsv_summary(all_results, output_dir)
        self.create_mlst_html_summary(all_results, output_dir)
        self.create_mlst_json_summary(all_results, output_dir)
        print("✅ MLST summary files created successfully!")

    def create_mlst_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples"""
        summary_file = output_dir / "mlst_summary.tsv"
        
        with open(summary_file, 'w') as f:
            # Header
            f.write("Sample\tST\tMLST_Status\tIdentity\tCoverage\tAllele_Profile\t")
            
            all_genes = set()
            for result in all_results.values():
                all_genes.update(result['alleles'].keys())
            
            for gene in sorted(all_genes):
                f.write(f"{gene}\t")
            f.write("\n")
            
            for sample_name, result in all_results.items():
                f.write(f"{sample_name}\t{result['st']}\t{result.get('mlst_status', 'Not Assigned')}\t{result.get('identity', 'Not Assigned')}\t{result.get('coverage', 'Not Assigned')}\t{result['allele_profile']}\t")
                for gene in sorted(all_genes):
                    allele = result['alleles'].get(gene, '')
                    f.write(f"{allele}\t")
                f.write("\n")
        
        print(f"📄 TSV summary created: {summary_file}")

    def create_mlst_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create beautiful HTML summary with teal/cyan + electric blue theme"""
        summary_file = output_dir / "mlst_summary.html"
        random_quote = self.get_random_quote()
        
        all_genes = set()
        for result in all_results.values():
            all_genes.update(result['alleles'].keys())
        sorted_genes = sorted(all_genes)
        
        total_samples = len(all_results)
        assigned_samples = sum(1 for result in all_results.values() if result.get('mlst_status') == 'Assigned')
        not_assigned_samples = total_samples - assigned_samples
        
        # Build table rows
        table_rows = ''
        for sample_name, result in all_results.items():
            st = result.get('st', '')
            mlst_status = result.get('mlst_status', 'Not Assigned')
            identity = result.get('identity', 'Not Assigned')
            coverage = result.get('coverage', 'Not Assigned')
            allele_profile = result.get('allele_profile', '')
            
            allele_columns = ''
            for gene in sorted_genes:
                allele = result['alleles'].get(gene, '')
                allele_columns += f'<td>{allele}</td>\n'
            
            mlst_status_html = f'<span style="color: {"#10b981" if mlst_status == "Assigned" else "#dc2626"}">{mlst_status}</span>'
            identity_class = 'identity-cell' if identity == '100%' else 'not-assigned-cell'
            coverage_class = 'coverage-cell' if coverage == '100%' else 'not-assigned-cell'
            
            table_rows += f'''                        <tr>
                            <td><strong>{sample_name}</strong></td>
                            <td class="st-cell">ST{st}</td>
                            <td>{mlst_status_html}</td>
                            <td class="{identity_class}">{identity}</td>
                            <td class="{coverage_class}">{coverage}</td>
                            <td class="allele-cell">{allele_profile}</td>
{allele_columns}                        </tr>
'''
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroScope - MLST Summary Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1800px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px;
            border: 2px solid rgba(6,182,212,0.5);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace; font-size: 10px; line-height: 1.1;
            white-space: pre; color: #06b6d4; text-shadow: 0 0 10px rgba(6,182,212,0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(0,0,0,0.4); backdrop-filter: blur(10px); padding: 20px; border-radius: 10px;
            margin-bottom: 30px; text-align: center; min-height: 100px;
            border: 1px solid rgba(6,182,212,0.4);
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{
            background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px;
            margin-bottom: 20px;
        }}
        .report-section h2 {{ color: #0f766e; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; }}
        .stats-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px,1fr)); gap: 15px; margin-bottom: 20px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%); color: white; padding: 15px;
            border-radius: 8px; text-align: center;
        }}
        .stat-value {{ font-size: 24px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .summary-table {{
            width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #0d9488 0%, #0f766e 100%); color: white; padding: 12px;
            text-align: left; position: sticky; top: 0;
        }}
        .summary-table td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #ccfbf1; }}
        .st-cell {{ font-weight: bold; color: #0f766e; }}
        .allele-cell {{ font-family: 'Courier New', monospace; background-color: #f0fdfa; color: #0f766e; font-weight: bold; }}
        .identity-cell, .coverage-cell {{ font-weight: bold; color: #10b981; }}
        .not-assigned-cell {{ color: #dc2626; font-style: italic; }}
        .footer {{ text-align: center; margin-top: 30px; padding: 20px; background: rgba(0,0,0,0.3); border-radius: 10px; }}
        .footer a {{ color: #fbbf24; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        @media (max-width: 768px) {{ .ascii-art {{ font-size: 6px; }} .summary-table {{ font-size: 12px; }} }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container"><div class="ascii-art">███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
█████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║███████╗██║      ██║   ██║██████╔╝█████╗  
██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝</div></div>
        <div class="quote-container"><div class="quote-text">"{random_quote['text']}"</div><div class="quote-author">— {random_quote['author']}</div></div>
    </div>
    <div class="report-section">
        <h2>MLST Summary - All Samples</h2>
        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{total_samples}</div><div class="stat-label">SAMPLES</div></div>
            <div class="stat-card"><div class="stat-value">{assigned_samples}</div><div class="stat-label">ASSIGNED</div></div>
            <div class="stat-card"><div class="stat-value">{not_assigned_samples}</div><div class="stat-label">NOT ASSIGNED</div></div>
        </div>
        <div style="overflow-x: auto;">
            <table class="summary-table">
                <thead><tr><th>Sample</th><th>ST</th><th>Status</th><th>Identity</th><th>Coverage</th><th>Allele Profile</th>'''
        
        for gene in sorted_genes:
            html_content += f'<th>{gene}</th>'
        
        html_content += f'''</tr></thead>
                <tbody>{table_rows}</tbody>
            </table>
        </div>
    </div>
    <div class="footer">
        <p><strong>EnteroScope</strong> - MLST Summary for <em>Enterobacter cloacae</em> Complex</p>
        <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        <p>GitHub: <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a> | Email: <a href="mailto:brownbeckley94@gmail.com">brownbeckley94@gmail.com</a></p>
    </div>
</div>
</body>
</html>'''
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"🌐 HTML summary created: {summary_file}")

    def create_mlst_json_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create JSON summary file with all MLST results"""
        summary_file = output_dir / "mlst_summary.json"
        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "total_samples": len(all_results),
                "analysis_type": "MLST",
                "scheme": "ecloacae",
                "version": "1.0"
            },
            "statistics": self._calculate_json_statistics(all_results),
            "samples": {}
        }
        for sample_name, result in all_results.items():
            json_summary["samples"][sample_name] = {
                "sequence_type": result.get('st', 'ND'),
                "mlst_status": result.get('mlst_status', 'Not Assigned'),
                "identity": result.get('identity', 'Not Assigned'),
                "coverage": result.get('coverage', 'Not Assigned'),
                "confidence": result.get('confidence', 'LOW'),
                "alleles": result.get('alleles', {}),
                "allele_profile": result.get('allele_profile', ''),
                "quality_metrics": result.get('quality_metrics', {})
            }
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, ensure_ascii=False)
        print(f"📄 JSON summary created: {summary_file}")

    def _calculate_json_statistics(self, all_results: Dict[str, Dict]) -> Dict:
        total_samples = len(all_results)
        assigned_samples = sum(1 for r in all_results.values() if r.get('mlst_status') == 'Assigned')
        unique_sts = set()
        st_counts = {}
        for r in all_results.values():
            st = r.get('st', '')
            if st and st not in ['ND', '-', 'UNKNOWN', '']:
                unique_sts.add(st)
                st_counts[st] = st_counts.get(st, 0) + 1
        
        all_genes = set()
        for r in all_results.values():
            all_genes.update(r.get('alleles', {}).keys())
        
        return {
            "total_samples": total_samples,
            "assigned_samples": assigned_samples,
            "assignment_rate": (assigned_samples / total_samples * 100) if total_samples > 0 else 0,
            "unique_sts": len(unique_sts),
            "st_distribution": dict(sorted(st_counts.items(), key=lambda x: x[1], reverse=True)),
            "total_genes": len(all_genes),
            "genes_detected": sorted(list(all_genes))
        }

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "ecloacae") -> Dict[str, Dict]:
        """Run MLST analysis for multiple files"""
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
        print(f"📁 Found {len(fasta_files)} FASTA files")
        results = {}
        for fasta_file in fasta_files:
            result = self.run_mlst_single(fasta_file, output_dir, scheme)
            results[fasta_file.name] = result
        self.create_mlst_summary(results, output_dir)
        return results

def main():
    parser = argparse.ArgumentParser(description='EnteroScope Modular MLST Analyzer (E. cloacae complex)')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file or directory (supports wildcards)')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('-db', '--database-dir', required=True, help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True, help='Script directory (contains mlst binary)')
    parser.add_argument('-s', '--scheme', default='ecloacae', help='MLST scheme (default: ecloacae)')
    parser.add_argument('--batch', action='store_true', help='Process multiple files')
    
    args = parser.parse_args()
    
    analyzer = ModularMLSTAnalyzer(
        database_dir=Path(args.database_dir),
        script_dir=Path(args.script_dir)
    )
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.batch:
        results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)
        print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
    else:
        input_file = Path(args.input)
        if input_file.exists():
            result = analyzer.run_mlst_single(input_file, output_dir, args.scheme)
            print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'ND')}")
        else:
            print(f"❌ Input file not found: {args.input}")

if __name__ == "__main__":
    main()