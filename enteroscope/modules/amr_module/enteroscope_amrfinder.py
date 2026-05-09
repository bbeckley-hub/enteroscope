#!/usr/bin/env python3
"""
EnteroScope AMRfinderPlus - BUNDLED VERSION with DYNAMIC DATABASE
Comprehensive AMR analysis for Enterobacter cloacae complex with HTML, TSV, and JSON reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
GitHub: https://github.com/bbeckley-hub/enteroscope
Date: 2026
Uses BUNDLED AMRFinderPlus 4.2.7 with LATEST DYNAMIC DATABASE
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import json
import random
from collections import defaultdict

class EnteroAMRfinderPlus:
    """AMRfinderPlus executor for Enterobacter cloacae complex with DYNAMIC database update"""
    
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        self.bundled_database = self._get_latest_database()
        
        if self.bundled_database is None:
            self.logger.warning("No AMRfinderPlus database found. Please run with --update-db to download.")
        
        db_version = self._get_database_version() if self.bundled_database else "Unknown"
        
        self.metadata = {
            "tool_name": "EnteroScope AMRfinderPlus",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub/enteroscope",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": db_version
        }
        
        # ========== Enterobacter cloacae complex specific gene sets ==========
        self.critical_carbapenemases = {
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaGES', 'blaIMI', 'blaSME', 'blaNMC', 'blaCcrA', 'blaCphA',
            'KPC', 'NDM', 'VIM', 'IMP', 'OXA-48', 'OXA-181', 'OXA-232',
            'GES', 'IMI', 'SME', 'NMC', 'CcrA', 'CphA'
        }
        
        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES', 'blaBEL',
            'CTX-M', 'SHV', 'TEM', 'PER', 'VEB', 'GES', 'BEL', 'BES',
            'blaCMY', 'blaDHA', 'blaFOX', 'blaMOX', 'blaACC', 'blaACT', 'blaMIR',
            'CMY', 'DHA', 'FOX', 'MOX', 'ACC', 'ACT', 'MIR'
        }
        
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'mcr1', 'mcr2', 'mcr3', 'mcr4', 'mcr5', 'mcr6', 'mcr7', 'mcr8', 'mcr9', 'mcr10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD', 'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }
        
        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH',
            'APH(3\')', 'APH(3\')-VI', 'APH(6)', 'AAC(3)', 'AAC(6\')', 'ANT(2")', 'ANT(3")', 'ANT(4")',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }
        
        self.high_risk_resistance = {
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA20', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
            'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
            'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
            'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19', 'qnrS1', 'qnrS2', 'qnrVC1', 'qnrVC4',
            'aac(6\')-Ib-cr', 'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
            'fosA', 'fosB', 'fosC', 'fosX',
            'arr-2', 'arr-3', 'arr-4', 'arr-5', 'arr-6', 'arr-7',
            'adeA', 'adeB', 'adeC', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'acrA', 'acrB', 'tolC', 'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF',
            'oprM', 'oprN', 'oprJ', 'mdtA', 'mdtB', 'mdtC', 'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY',
            'mdx', 'acrAB', 'acrEF', 'mdfA', 'cmlA'
        }
        
        self.high_risk_genes = (self.critical_carbapenemases.union(self.critical_esbls)
                                .union(self.critical_colistin).union(self.critical_aminoglycoside)
                                .union(self.high_risk_resistance))
        self.critical_risk_genes = self.critical_carbapenemases.union(self.critical_colistin)
        
        # Science quotes
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
            {"text": "EnteroScope makes advanced genomic typing for Enterobacter cloacae complex accessible to all.", "author": "Brown Beckley"}
        ]
        
        # EnteroScope ASCII art (teal/cyan coloured in CSS)
        self.ascii_art = r"""
███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
█████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║███████╗██║      ██║   ██║██████╔╝█████╗  
██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝
"""
    
    def get_random_quote(self):
        return random.choice(self.science_quotes)
    
    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> int:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except:
            return 8
    
    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            user_cpus = min(user_cpus, 64)
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            total = min(total, 64)
            if total <= 4:
                opt = total
            elif total <= 8:
                opt = total - 1
            elif total <= 16:
                opt = max(8, total - 2)
            elif total <= 32:
                opt = max(16, total - 3)
            else:
                opt = min(64, int(total * 0.95))
            opt = max(1, min(opt, total, 64))
            self._log_resource_info(opt, total)
            return opt
        except:
            return min(os.cpu_count() or 4, 64)
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            self.logger.info(f"Using CPU cores: {cpus} ({cpus/total_cores*100:.1f}%)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus >= 48:
            self.logger.info("💡 Performance: EXTREME SPEED MODE 🚀")
        elif cpus >= 32:
            self.logger.info("💡 Performance: ULTRA MAXIMUM SPEED MODE 🚀")
        elif cpus >= 16:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")
        elif cpus >= 8:
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: Standard multi-core")
        self.logger.info("📝 STRATEGY: Processing multiple Enterobacter genomes concurrently with optimal core allocation")
    
    def _get_latest_database(self) -> Optional[str]:
        db_root = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.exists(db_root):
            return None
        candidates = []
        for item in os.listdir(db_root):
            full = os.path.join(db_root, item)
            if os.path.isdir(full) and item.startswith('20'):
                candidates.append(item)
        if not candidates:
            return None
        latest = sorted(candidates)[-1]
        return os.path.join(db_root, latest)
    
    def _get_database_version(self) -> str:
        if not self.bundled_database:
            return "Unknown"
        version_file = os.path.join(self.bundled_database, "version.txt")
        if os.path.exists(version_file):
            with open(version_file, 'r') as f:
                return f.read().strip()
        return os.path.basename(self.bundled_database)
    
    def update_database(self, force: bool = False) -> bool:
        if not os.path.exists(self.bundled_update):
            self.logger.error(f"amrfinder_update not found at {self.bundled_update}")
            return False
        if not os.access(self.bundled_update, os.X_OK):
            os.chmod(self.bundled_update, 0o755)
        db_dir = os.path.join(self.module_dir, "data", "amrfinder_db")
        os.makedirs(db_dir, exist_ok=True)
        cmd = [self.bundled_update, "--database", db_dir]
        if force:
            cmd.append("--force_update")
            self.logger.info("🔨 Force update mode: existing database will be overwritten")
        self.logger.info("Updating AMRfinderPlus database...")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info("Database update completed successfully.")
            self.bundled_database = self._get_latest_database()
            if self.bundled_database:
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info(f"New database version: {self.metadata['database_version']}")
                return True
            else:
                self.logger.error("Database update succeeded but no database folder found.")
                return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Database update failed: {e.stderr}")
            return False
    
    def check_amrfinder_installed(self) -> bool:
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
                return False
            if not os.access(self.bundled_amrfinder, os.X_OK):
                os.chmod(self.bundled_amrfinder, 0o755)
            result = subprocess.run([self.bundled_amrfinder, '--version'], capture_output=True, text=True, check=True)
            self.logger.info(f"Bundled AMRfinderPlus version: {result.stdout.strip()}")
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info(f"✅ Bundled database found: {self.bundled_database}")
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info(f"✅ Database version: {self.metadata['database_version']}")
                return True
            else:
                self.logger.warning("⚠️ Bundled database not found. Please run --update-db to download.")
                return False
        except Exception as e:
            self.logger.error(f"Bundled AMRfinderPlus check failed: {e}")
            return False
    
    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        run_cpus = self.cpus
        cmd = [
            self.bundled_amrfinder,
            '-n', genome_file,
            '--output', output_file,
            '--threads', str(run_cpus),
            '--plus',
            '--organism', 'Enterobacter_cloacae'
        ]
        if self.bundled_database and os.path.exists(self.bundled_database):
            cmd.extend(['--database', self.bundled_database])
            self.logger.info(f"Using bundled database: {self.bundled_database}")
        else:
            self.logger.warning("Using default AMRfinderPlus database location")
        self.logger.info(f"Running AMRfinderPlus: {genome_name} (using {run_cpus} CPU cores)")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = self._parse_amrfinder_output(output_file)
            self._create_amrfinder_html_report(genome_name, hits, output_dir)
            return {'genome': genome_name, 'output_file': output_file, 'hits': hits, 'hit_count': len(hits), 'status': 'success'}
        except subprocess.CalledProcessError as e:
            self.logger.error(f"AMRfinderPlus failed for {genome_name}: {e.stderr}")
            return {'genome': genome_name, 'output_file': output_file, 'hits': [], 'hit_count': 0, 'status': 'failed', 'error': e.stderr}
    
    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
            if not lines or len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = dict(zip(headers, parts))
                    processed = {
                        'protein_id': hit.get('Protein id', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': hit.get('Element symbol', ''),
                        'sequence_name': hit.get('Element name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Type', ''),
                        'element_subtype': hit.get('Subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                        'identity': hit.get('% Identity to reference', '').replace('%', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Closest reference accession', ''),
                        'closest_name': hit.get('Closest reference name', ''),
                        'hmm_id': hit.get('HMM accession', ''),
                        'hmm_description': hit.get('HMM description', '')
                    }
                    hits.append(processed)
                else:
                    self.logger.warning(f"Line {line_num} has {len(parts)} parts, expected {len(headers)}: {line[:100]}")
        except Exception as e:
            self.logger.error(f"Error parsing {amrfinder_file}: {e}")
        self.logger.info(f"Parsed {len(hits)} AMR hits from {amrfinder_file}")
        return hits
    
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        analysis = self._analyze_enterobacter_amr_results(hits)
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Teal/Cyan CSS
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroScope - AMRfinderPlus Analysis Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7); padding:20px; border-radius:15px; margin-bottom:20px;
            border: 2px solid #06b6d4;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace; font-size:10px; white-space:pre;
            color: #06b6d4; text-shadow: 0 0 10px #06b6d4; overflow-x:auto;
        }}
        .quote-container {{
            background: rgba(0,0,0,0.4); backdrop-filter:blur(10px); padding:20px; border-radius:10px;
            margin-bottom:30px; text-align:center; border-left:4px solid #06b6d4;
            transition:opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; color:#fff; }}
        .quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
        .report-section {{ background: rgba(255,255,255,0.95); color:#1f2937; padding:25px; border-radius:12px; margin-bottom:20px; box-shadow:0 4px 15px rgba(0,0,0,0.2); }}
        .report-section h2 {{ color:#0f766e; border-bottom:3px solid #3b82f6; padding-bottom:10px; margin-bottom:20px; font-size:24px; }}
        .metrics-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(250px,1fr)); gap:20px; margin-top:15px; }}
        .metric-card {{ background:linear-gradient(135deg,#06b6d4 0%,#0891b2 100%); color:white; padding:20px; border-radius:8px; box-shadow:0 4px 12px rgba(0,0,0,0.15); }}
        .metric-label {{ font-size:14px; opacity:0.9; margin-bottom:5px; }}
        .metric-value {{ font-size:24px; font-weight:bold; }}
        .summary-table {{ width:100%; border-collapse:collapse; margin-top:20px; font-size:14px; }}
        .summary-table th {{ background:linear-gradient(135deg,#3b82f6 0%,#2563eb 100%); color:white; padding:12px; text-align:left; }}
        .summary-table td {{ padding:12px; border-bottom:1px solid #e5e7eb; }}
        .summary-table tr:nth-child(even) {{ background-color:#f8fafc; }}
        .summary-table tr:hover {{ background-color:#e0f2fe; }}
        .table-responsive {{ overflow-x:auto; margin:20px 0; }}
        .present {{ background-color:#d4edda; }}
        .high-risk {{ background-color:#fef3c7; }}
        .critical {{ background-color:#fee2e2; font-weight:bold; }}
        .footer {{ text-align:center; margin-top:30px; padding:20px; background:rgba(0,0,0,0.3); border-radius:10px; font-size:14px; }}
        .footer a {{ color:#fbbf24; text-decoration:none; }}
        .footer a:hover {{ text-decoration:underline; }}
        .timestamp {{ color:#fbbf24; font-weight:bold; }}
        .authorship {{ margin-top:15px; padding:15px; background:rgba(255,255,255,0.1); border-radius:8px; font-size:12px; }}
        .risk-badge {{ display:inline-block; background:#dc2626; color:white; padding:5px 10px; border-radius:15px; margin:2px; font-size:0.9em; }}
        .warning-badge {{ background:#f59e0b; color:black; }}
        .product-cell {{ white-space:normal !important; word-wrap:break-word; max-width:500px; min-width:200px; }}
        @media (max-width:768px) {{ .ascii-art {{ font-size:6px; }} .metrics-grid {{ grid-template-columns:1fr; }} .summary-table {{ font-size:12px; }} }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
        <div class="quote-container" id="quoteContainer">
            <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
            <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
        </div>
    </div>
    <div class="report-section">
        <h2>📊 Enterobacter cloacae AMR Summary</h2>
        <div class="metrics-grid">
            <div class="metric-card"><div class="metric-label">Total AMR Genes</div><div class="metric-value">{analysis['total_genes']}</div></div>
            <div class="metric-card"><div class="metric-label">High Risk Genes</div><div class="metric-value">{analysis['high_risk_genes']}</div></div>
            <div class="metric-card" style="background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%);"><div class="metric-label">Critical Risk</div><div class="metric-value">{analysis['critical_risk_genes']}</div></div>
            <div class="metric-card"><div class="metric-label">Analysis Date</div><div class="metric-value">{datetime.now().strftime('%Y-%m-%d')}</div></div>
        </div>
        <p><strong>Genome:</strong> {genome_name}</p>
        <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
        <p><strong>Database:</strong> {self.metadata['database_version']}</p>
    </div>
"""
        if analysis['critical_risk_genes'] > 0:
            html_content += f"""
    <div class="report-section" style="border-left:4px solid #dc2626;">
        <h2 style="color:#dc2626;">⚠️ CRITICAL RISK AMR GENES DETECTED</h2>
        <p><strong>{analysis['critical_risk_genes']} critical risk antimicrobial resistance genes found:</strong></p>
        <div style="margin:10px 0;">"""
            for gene in analysis['critical_risk_list']:
                html_content += f'<span class="risk-badge">🚨 {gene}</span>'
            html_content += "</div></div>"
        
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            html_content += f"""
    <div class="report-section" style="border-left:4px solid #f59e0b;">
        <h2 style="color:#f59e0b;">⚠️ High-Risk AMR Genes Detected</h2>
        <p><strong>{analysis['high_risk_genes']} high-risk antimicrobial resistance genes found:</strong></p>
        <div style="margin:10px 0;">"""
            for gene in analysis['high_risk_list']:
                html_content += f'<span class="warning-badge">{gene}</span>'
            html_content += "</div></div>"
        
        if any(analysis['resistance_mechanisms'].values()):
            html_content += """
    <div class="report-section">
        <h2 style="color:#0f766e; border-bottom:3px solid #3b82f6;">🔬 Resistance Mechanism Breakdown</h2>"""
            mech = analysis['resistance_mechanisms']
            if mech['carbapenemase']:
                html_content += f'<div style="margin:10px 0; padding:10px; background:#fee2e2; border-radius:5px;"><strong>Carbapenemase Genes (CRITICAL):</strong> {", ".join(mech["carbapenemase"])}</div>'
            if mech['esbl']:
                html_content += f'<div style="margin:10px 0; padding:10px; background:#fff3cd; border-radius:5px;"><strong>ESBL Genes:</strong> {", ".join(mech["esbl"])}</div>'
            if mech['colistin_resistance']:
                html_content += f'<div style="margin:10px 0; padding:10px; background:#fee2e2; border-radius:5px;"><strong>Colistin Resistance (CRITICAL):</strong> {", ".join(mech["colistin_resistance"])}</div>'
            if mech['aminoglycoside_resistance']:
                html_content += f'<div style="margin:10px 0; padding:10px; background:#d1ecf1; border-radius:5px;"><strong>Aminoglycoside Resistance:</strong> {", ".join(mech["aminoglycoside_resistance"])}</div>'
            if mech['efflux_pumps']:
                html_content += f'<div style="margin:10px 0; padding:10px; background:#e2e3e5; border-radius:5px;"><strong>Efflux Pumps:</strong> {", ".join(mech["efflux_pumps"])}</div>'
            if mech['other_amr']:
                html_content += f'<div style="margin:10px 0; padding:10px; background:#f8f9fa; border-radius:5px;"><strong>Other AMR Genes:</strong> {", ".join(mech["other_amr"])}</div>'
            html_content += "</div>"
        
        if analysis['resistance_classes']:
            html_content += """
    <div class="report-section">
        <h2 style="color:#0f766e; border-bottom:3px solid #3b82f6;">🧪 Resistance Classes Detected</h2>
        <div class="table-responsive"><table class="summary-table"><thead><tr><th>Resistance Class</th><th>Gene Count</th><th>Genes</th></tr></thead><tbody>"""
            for cls, genes in analysis['resistance_classes'].items():
                html_content += f"<tr><td><strong>{cls}</strong></td><td>{len(genes)}</td><td class='product-cell'>{', '.join(genes)}</td></tr>"
            html_content += "</tbody></table></div></div>"
        
        if hits:
            html_content += """
    <div class="report-section">
        <h2 style="color:#0f766e; border-bottom:3px solid #3b82f6;">🔬 Detailed AMR Genes Detected</h2>
        <div class="table-responsive"><table class="summary-table"><thead><tr><th>Gene Symbol</th><th>Sequence Name</th><th>Class</th><th>Subclass</th><th>Coverage</th><th>Identity</th><th>Scope</th></tr></thead><tbody>"""
            for hit in hits:
                gene = hit.get('gene_symbol', '')
                row_class = "critical" if gene in analysis['critical_risk_list'] else ("high-risk" if gene in analysis['high_risk_list'] else "present")
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{gene}</strong></td>
                        <td class="product-cell">{hit.get('sequence_name', '')}</td>
                        <td>{hit.get('class', '')}</td>
                        <td>{hit.get('subclass', '')}</td>
                        <td>{hit.get('coverage', '')}%</td>
                        <td>{hit.get('identity', '')}%</td>
                        <td>{hit.get('scope', '')}</td>
                    </tr>"""
            html_content += "</tbody></table></div></div>"
        else:
            html_content += '<div class="report-section"><h2 style="color:#0f766e;">✅ No AMR Genes Detected</h2><p>No antimicrobial resistance genes found in this E. cloacae genome.</p></div>'
        
        html_content += f"""
    <div class="footer">
        <p><strong>EnteroScope</strong> - AMRfinderPlus Analysis Module</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p><strong>Technical Support & Inquiries:</strong></p>
            <p>Author: Brown Beckley | GitHub: <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a></p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
<script>
    const quotes = {json.dumps([q['text'] for q in self.science_quotes])};
    const authors = {json.dumps([q['author'] for q in self.science_quotes])};
    let idx = 0;
    function rotate() {{
        const container = document.getElementById('quoteContainer');
        const textDiv = document.getElementById('quoteText');
        const authorDiv = document.getElementById('quoteAuthor');
        container.style.opacity = '0';
        setTimeout(() => {{
            textDiv.innerHTML = '"' + quotes[idx] + '"';
            authorDiv.innerHTML = '— ' + authors[idx];
            container.style.opacity = '1';
            idx = (idx + 1) % quotes.length;
        }}, 500);
    }}
    setInterval(rotate, 10000);
    document.addEventListener('DOMContentLoaded', rotate);
</script>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info(f"EnteroScope AMRfinderPlus HTML report generated: {html_file}")
    
    def _analyze_enterobacter_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'carbapenemase': [], 'esbl': [], 'colistin_resistance': [],
                'aminoglycoside_resistance': [], 'efflux_pumps': [], 'other_amr': []
            }
        }
        for hit in hits:
            gene = hit.get('gene_symbol', '')
            rclass = hit.get('class', '')
            if not gene:
                continue
            self._categorize_enterobacter_mechanism(gene, analysis)
            if gene in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene)
            if gene in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene)
            if rclass:
                analysis['resistance_classes'].setdefault(rclass, [])
                if gene not in analysis['resistance_classes'][rclass]:
                    analysis['resistance_classes'][rclass].append(gene)
        return analysis
    
    def _categorize_enterobacter_mechanism(self, gene: str, analysis: Dict):
        gl = gene.lower()
        if any(c in gl for c in ['kpc', 'ndm', 'vim', 'imp', 'oxa-48', 'ges', 'imi', 'sme', 'nmc']):
            if any(c in gl for c in [x.lower() for x in self.critical_carbapenemases]):
                analysis['resistance_mechanisms']['carbapenemase'].append(gene)
                return
        if any(e in gl for e in ['ctx-m', 'shv', 'tem', 'per', 'veb', 'ges', 'cmy', 'dha', 'fox']):
            if any(c in gl for c in [x.lower() for x in self.critical_esbls]):
                analysis['resistance_mechanisms']['esbl'].append(gene)
                return
        if 'mcr' in gl:
            if any(c in gl for c in [x.lower() for x in self.critical_colistin]):
                analysis['resistance_mechanisms']['colistin_resistance'].append(gene)
                return
        if any(ag in gl for ag in ['arm', 'rmt', 'aph', 'aac', 'aad', 'str', 'ant']):
            if any(c in gl for c in [x.lower() for x in self.critical_aminoglycoside]):
                analysis['resistance_mechanisms']['aminoglycoside_resistance'].append(gene)
                return
        if any(ef in gl for ef in ['ade', 'abe', 'mex', 'acr', 'emr', 'mdt', 'tolc']):
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene)
            return
        analysis['resistance_mechanisms']['other_amr'].append(gene)
    
    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating EnteroScope AMR summary files...")
        summary_file = os.path.join(output_base, "enteroscope_amrfinder_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("Genome\tGene_Symbol\tSequence_Name\tClass\tSubclass\tCoverage\tIdentity\tScope\tElement_Type\tAccession\tContig\tStart\tStop\n")
            for gname, res in all_results.items():
                for hit in res['hits']:
                    row = [gname, hit.get('gene_symbol',''), hit.get('sequence_name',''), hit.get('class',''), hit.get('subclass',''),
                           hit.get('coverage',''), hit.get('identity',''), hit.get('scope',''), hit.get('element_type',''),
                           hit.get('accession',''), hit.get('contig_id',''), hit.get('start',''), hit.get('stop','')]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"✓ TSV summary: {summary_file}")
        stats_file = os.path.join(output_base, "enteroscope_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\tResistance_Classes\tGene_List\n")
            for gname, res in all_results.items():
                genes = list(set(h.get('gene_symbol','') for h in res['hits'] if h.get('gene_symbol')))
                high = sum(1 for g in genes if g in self.high_risk_genes)
                crit = sum(1 for g in genes if g in self.critical_risk_genes)
                classes = list(set(h.get('class','') for h in res['hits'] if h.get('class')))
                f.write(f"{gname}\t{res['hit_count']}\t{high}\t{crit}\t{','.join(classes)}\t{','.join(genes)}\n")
        self.logger.info(f"✓ Statistics summary: {stats_file}")
        self.create_amr_json_summaries(all_results, output_base)
        self.create_amr_master_json_summary(all_results, output_base)
        self._create_summary_html_report(all_results, output_base)
    
    def create_amr_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating JSON AMR summaries...")
        all_hits = []
        for gname, res in all_results.items():
            for hit in res['hits']:
                h = hit.copy()
                h['genome'] = gname
                all_hits.append(h)
        if not all_hits:
            self.logger.info("No AMR hits found, skipping JSON summaries")
            return
        gene_freq = {}
        for hit in all_hits:
            gene = hit.get('gene_symbol')
            if not gene:
                continue
            if gene not in gene_freq:
                gene_freq[gene] = {'count':0, 'genomes':set(), 'details':[]}
            gene_freq[gene]['count'] += 1
            gene_freq[gene]['genomes'].add(hit['genome'])
            gene_freq[gene]['details'].append({'genome':hit['genome'], 'sequence_name':hit.get('sequence_name',''), 'class':hit.get('class',''), 'subclass':hit.get('subclass',''), 'coverage':hit.get('coverage',''), 'identity':hit.get('identity',''), 'accession':hit.get('accession','')})
        for g in gene_freq:
            gene_freq[g]['genomes'] = list(gene_freq[g]['genomes'])
        json_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_hits': len(all_hits),
                'total_genomes': len(all_results),
                'unique_genes': len(gene_freq),
                'critical_risk_genes_found': list(self.critical_risk_genes.intersection(gene_freq.keys())),
                'high_risk_genes_found': list(self.high_risk_genes.intersection(gene_freq.keys()))
            },
            'gene_frequency': gene_freq,
            'summary_by_genome': self._create_amr_genome_summary(all_results),
            'hits': all_hits[:1000]
        }
        json_file = os.path.join(output_base, "enteroscope_amrfinder_summary.json")
        with open(json_file, 'w') as f:
            json.dump(json_summary, f, indent=2)
        self.logger.info(f"✓ JSON summary: {json_file}")
    
    def _create_amr_genome_summary(self, all_results: Dict) -> Dict:
        summary = {}
        for gname, res in all_results.items():
            summary[gname] = {'total_genes': res['hit_count'], 'genes': set(), 'critical_risk_genes': [], 'high_risk_genes': [], 'resistance_classes': set()}
            for hit in res['hits']:
                gene = hit.get('gene_symbol')
                if not gene:
                    continue
                summary[gname]['genes'].add(gene)
                if hit.get('class'):
                    summary[gname]['resistance_classes'].add(hit['class'])
                if gene in self.critical_risk_genes:
                    if gene not in summary[gname]['critical_risk_genes']:
                        summary[gname]['critical_risk_genes'].append(gene)
                elif gene in self.high_risk_genes:
                    if gene not in summary[gname]['high_risk_genes']:
                        summary[gname]['high_risk_genes'].append(gene)
        for g in summary:
            summary[g]['genes'] = list(summary[g]['genes'])
            summary[g]['resistance_classes'] = list(summary[g]['resistance_classes'])
        return summary
    
    def create_amr_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating master JSON summary...")
        master = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
                'critical_risk_genes': list(self.critical_risk_genes),
                'high_risk_genes': list(self.high_risk_genes)
            },
            'genome_summaries': {},
            'critical_findings': {},
            'cross_genome_patterns': {}
        }
        genomes_crit = 0
        genomes_high = 0
        all_hits_by_gene = {}
        for gname, res in all_results.items():
            genes = [h.get('gene_symbol') for h in res['hits'] if h.get('gene_symbol')]
            unique = set(genes)
            classes = {h.get('class') for h in res['hits'] if h.get('class')}
            has_crit = any(g in self.critical_risk_genes for g in unique)
            has_high = any(g in self.high_risk_genes for g in unique)
            if has_crit:
                genomes_crit += 1
            if has_high:
                genomes_high += 1
            master['genome_summaries'][gname] = {
                'total_hits': res['hit_count'],
                'unique_genes': len(unique),
                'has_critical_risk': has_crit,
                'has_high_risk': has_high,
                'resistance_classes': list(classes),
                'genes': list(unique)
            }
            if has_crit:
                master['critical_findings'].setdefault('genomes_with_critical', []).append({'genome': gname, 'critical_genes': [g for g in unique if g in self.critical_risk_genes]})
            for g in unique:
                if g not in all_hits_by_gene:
                    all_hits_by_gene[g] = {'count':0, 'genomes':[], 'classes':set()}
                all_hits_by_gene[g]['count'] += 1
                if gname not in all_hits_by_gene[g]['genomes']:
                    all_hits_by_gene[g]['genomes'].append(gname)
                if classes:
                    all_hits_by_gene[g]['classes'].update(classes)
        for g in all_hits_by_gene:
            all_hits_by_gene[g]['classes'] = list(all_hits_by_gene[g]['classes'])
        common = {g: data for g, data in all_hits_by_gene.items() if data['count'] > 1}
        top = sorted([(g, data) for g, data in all_hits_by_gene.items()], key=lambda x: x[1]['count'], reverse=True)[:20]
        master['cross_genome_patterns'] = {
            'total_genes_found': len(all_hits_by_gene),
            'genomes_with_critical_risk': genomes_crit,
            'genomes_with_high_risk': genomes_high,
            'common_genes': common,
            'top_genes': top,
            'risk_gene_distribution': {
                'critical_risk_found': len([g for g in all_hits_by_gene if g in self.critical_risk_genes]),
                'high_risk_found': len([g for g in all_hits_by_gene if g in self.high_risk_genes]),
                'standard_risk_found': len([g for g in all_hits_by_gene if g not in self.high_risk_genes])
            }
        }
        json_file = os.path.join(output_base, "enteroscope_amrfinder_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master, f, indent=2)
        self.logger.info(f"✓ Master JSON: {json_file}")
    
    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        # Collect data
        total_genomes = len(all_results)
        total_hits = sum(r['hit_count'] for r in all_results.values())
        genes_per_genome = {}
        gene_freq = defaultdict(set)
        critical_genes = set()
        high_risk_genes = set()
        genomes_with_critical = 0
        genomes_with_high = 0
        for gname, res in all_results.items():
            genes = set()
            for hit in res['hits']:
                g = hit.get('gene_symbol')
                if g:
                    genes.add(g)
                    gene_freq[g].add(gname)
                    if g in self.critical_risk_genes:
                        critical_genes.add(g)
                    elif g in self.high_risk_genes:
                        high_risk_genes.add(g)
            genes_per_genome[gname] = genes
            if any(g in self.critical_risk_genes for g in genes):
                genomes_with_critical += 1
            if any(g in self.high_risk_genes for g in genes):
                genomes_with_high += 1
        
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Build HTML (teal/cyan theme)
        html_parts = []
        html_parts.append(f"""<!DOCTYPE html>
<html>
<head><title>EnteroScope - AMRfinderPlus Batch Summary</title>
<meta charset="UTF-8">
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color: #ffffff;
    padding: 20px;
    min-height: 100vh;
}}
.container {{ max-width: 1400px; margin: 0 auto; }}
.header {{ text-align: center; margin-bottom: 30px; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px; margin-bottom:20px;
    border: 2px solid #06b6d4;
}}
.ascii-art {{
    font-family: 'Courier New', monospace; font-size:10px; white-space:pre;
    color: #06b6d4; text-shadow: 0 0 10px #06b6d4; overflow-x:auto;
}}
.quote-container {{
    background: rgba(0,0,0,0.4); backdrop-filter:blur(10px); padding:20px; border-radius:10px;
    margin-bottom:30px; text-align:center; border-left:4px solid #06b6d4;
}}
.quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; }}
.quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
.report-section {{
    background: rgba(255,255,255,0.95); color:#1f2937; padding:25px; border-radius:12px;
    margin-bottom:20px; box-shadow:0 4px 15px rgba(0,0,0,0.2);
}}
.report-section h2 {{ color:#0f766e; border-bottom:3px solid #3b82f6; padding-bottom:10px; margin-bottom:20px; }}
.metrics-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr)); gap:20px; margin:20px 0; }}
.metric-card {{ background:linear-gradient(135deg,#06b6d4 0%,#0891b2 100%); color:white; padding:20px; border-radius:8px; text-align:center; }}
.summary-table {{ width:100%; border-collapse:collapse; margin-top:20px; font-size:14px; }}
.summary-table th {{ background:linear-gradient(135deg,#3b82f6 0%,#2563eb 100%); color:white; padding:12px; text-align:left; }}
.summary-table td {{ padding:12px; border-bottom:1px solid #e5e7eb; }}
.summary-table tr:nth-child(even) {{ background-color:#f8fafc; }}
.summary-table tr:hover {{ background-color:#e0f2fe; }}
.risk-badge {{ display:inline-block; background:#dc2626; color:white; padding:5px 10px; border-radius:15px; margin:2px; }}
.warning-badge {{ background:#f59e0b; color:black; }}
.footer {{ text-align:center; margin-top:30px; padding:20px; background:rgba(0,0,0,0.3); border-radius:10px; }}
.footer a {{ color:#fbbf24; text-decoration:none; }}
.timestamp {{ color:#fbbf24; }}
.sequence-cell {{ white-space:normal; word-wrap:break-word; max-width:300px; }}
</style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
        <div class="quote-container">
            <div class="quote-text" id="summary-quote">"{random_quote['text']}"</div>
            <div class="quote-author">— {random_quote['author']}</div>
        </div>
    </div>
    <div class="report-section">
        <h2>📊 Batch AMRfinderPlus Summary</h2>
        <div class="metrics-grid">
            <div class="metric-card"><div class="metric-value">{total_genomes}</div><div>Total Genomes</div></div>
            <div class="metric-card"><div class="metric-value">{total_hits}</div><div>Total AMR Genes</div></div>
            <div class="metric-card" style="background:#dc2626;"><div class="metric-value">{genomes_with_critical}</div><div>Critical Risk Genomes</div></div>
        </div>
        <p>Tool: {self.metadata['version']} | AMRfinderPlus: {self.metadata['amrfinder_version']} | Database: {self.metadata['database_version']}</p>
    </div>
""")
        
        if critical_genes:
            html_parts.append(f'<div class="report-section" style="border-left:4px solid #dc2626;"><h2 style="color:#dc2626;">🚨 CRITICAL RISK AMR GENES ACROSS ALL GENOMES</h2><p>{len(critical_genes)} unique critical risk genes found in {genomes_with_critical} genomes:</p><div>')
            for g in sorted(critical_genes):
                html_parts.append(f'<span class="risk-badge">🚨 {g}</span>')
            html_parts.append('</div></div>')
        elif high_risk_genes:
            html_parts.append(f'<div class="report-section" style="border-left:4px solid #f59e0b;"><h2 style="color:#f59e0b;">⚠️ High-Risk AMR Genes Detected</h2><p>{len(high_risk_genes)} unique high-risk genes found across {genomes_with_high} genomes:</p><div>')
            for g in sorted(high_risk_genes):
                html_parts.append(f'<span class="warning-badge">{g}</span>')
            html_parts.append('</div></div>')
        
        html_parts.append('<div class="report-section"><h2>🔍 Genes by Genome</h2><div class="table-responsive"><table class="summary-table"><thead><tr><th>Genome</th><th>Gene Count</th><th>Critical Genes</th><th>High Risk Genes</th></tr></thead><tbody>')
        for genome in sorted(genes_per_genome):
            genes = genes_per_genome[genome]
            crit = [g for g in genes if g in self.critical_risk_genes]
            high = [g for g in genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            html_parts.append(f'<tr><td><strong>{genome}</strong></td><td>{len(genes)}</td><td class="sequence-cell">{", ".join(crit) if crit else "None"}</td><td class="sequence-cell">{", ".join(high) if high else "None"}</td></tr>')
        html_parts.append('</tbody></table></div></div>')
        
        html_parts.append('<div class="report-section"><h2>📈 Gene Frequency</h2><div class="table-responsive"><table class="summary-table"><thead><tr><th>Gene</th><th>Frequency</th><th>Prevalence</th><th>Risk Level</th><th>Genomes</th></tr></thead><tbody>')
        for gene, genomes in sorted(gene_freq.items(), key=lambda x: len(x[1]), reverse=True):
            freq = len(genomes)
            percent = freq / total_genomes * 100
            if gene in self.critical_risk_genes:
                risk = '<span class="risk-badge">CRITICAL</span>'
            elif gene in self.high_risk_genes:
                risk = '<span class="warning-badge">HIGH</span>'
            else:
                risk = '<span class="success-badge">Standard</span>'
            if percent >= 75:
                badge = '<span class="risk-badge">Very High</span>'
            elif percent >= 50:
                badge = '<span class="warning-badge">High</span>'
            elif percent >= 25:
                badge = '<span class="success-badge">Medium</span>'
            elif percent >= 10:
                badge = '<span class="success-badge">Low</span>'
            else:
                badge = '<span class="success-badge">Rare</span>'
            html_parts.append(f'<tr><td><strong>{gene}</strong></td><td>{freq} ({percent:.1f}%)</td><td>{badge}</td><td>{risk}</td><td class="sequence-cell">{", ".join(sorted(genomes))}</td></tr>')
        html_parts.append('</tbody></table></div></div>')
        
        html_parts.append(f"""
    <div class="footer">
        <p><strong>EnteroScope</strong> - AMRfinderPlus Batch Summary</p>
        <p class="timestamp">Generated: {current_time}</p>
        <p><a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a> | brownbeckley94@gmail.com</p>
    </div>
</div>
<script>
    const quotes = {json.dumps([q['text'] for q in self.science_quotes])};
    const authors = {json.dumps([q['author'] for q in self.science_quotes])};
    let idx = 0;
    function rotate() {{
        const qdiv = document.getElementById('summary-quote');
        const adiv = document.querySelector('.quote-author');
        if (qdiv && adiv) {{
            qdiv.innerHTML = '"' + quotes[idx] + '"';
            adiv.innerHTML = '— ' + authors[idx];
            idx = (idx + 1) % quotes.length;
        }}
    }}
    setInterval(rotate, 10000);
    document.addEventListener('DOMContentLoaded', rotate);
</script>
</body>
</html>""")
        
        html_file = os.path.join(output_base, "enteroscope_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write('\n'.join(html_parts))
        self.logger.info(f"✓ Summary HTML report: {html_file}")
    
    def process_single_genome(self, genome_file: str, output_base: str = "enteroscope_amrfinder_results") -> Dict:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        result = self.run_amrfinder_single_genome(genome_file, results_dir)
        self.logger.info(f"{'✓' if result['status']=='success' else '✗'} {genome_name}: {result['hit_count']} AMR hits")
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "enteroscope_amrfinder_results") -> Dict:
        if not self.check_amrfinder_installed():
            raise RuntimeError("BUNDLED AMRfinderPlus not properly installed")
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa", f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        files = []
        for p in fasta_patterns:
            files.extend(glob.glob(p))
        files = list(set(files))
        if not files:
            raise FileNotFoundError(f"No FASTA files matching '{genome_pattern}'")
        self.logger.info(f"Found {len(files)} genomes: {[Path(f).name for f in files]}")
        os.makedirs(output_base, exist_ok=True)
        max_concurrent = max(1, min(self.cpus, len(files), int(self.available_ram / 1.5)))
        self.logger.info(f"🚀 MAXIMUM SPEED: Using {max_concurrent} concurrent jobs, each using {self.cpus} threads")
        all_results = {}
        with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
            futures = {ex.submit(self.process_single_genome, f, output_base): f for f in files}
            for fut in as_completed(futures):
                try:
                    res = fut.result()
                    all_results[res['genome']] = res
                except Exception as e:
                    self.logger.error(f"Failed {futures[fut]}: {e}")
        self.create_amr_summary(all_results, output_base)
        self.logger.info(f"=== ANALYSIS COMPLETE: {len(all_results)} genomes processed ===")
        return all_results

def main():
    parser = argparse.ArgumentParser(description='EnteroScope AMRfinderPlus - Enterobacter cloacae complex AMR analysis')
    parser.add_argument('pattern', nargs='?', help='File pattern for genomes (e.g., "*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, help='Number of CPU cores (max 64)')
    parser.add_argument('--output', '-o', default='enteroscope_amrfinder_results', help='Output directory')
    parser.add_argument('--update-db', action='store_true', help='Update AMRfinderPlus database to latest version and exit')
    parser.add_argument('--force-update', action='store_true', help='Force update (overwrite existing database)')
    parser.add_argument('--db-version', action='store_true', help='Show current database version and exit')
    args = parser.parse_args()
    
    if args.update_db or args.force_update or args.db_version:
        executor = EnteroAMRfinderPlus(cpus=args.cpus)
        if args.update_db or args.force_update:
            print("Updating AMRfinderPlus database...")
            success = executor.update_database(force=args.force_update)
            sys.exit(0 if success else 1)
        if args.db_version:
            print(f"Database version: {executor.metadata['database_version']}")
            print(f"Database path: {executor.bundled_database or 'Not found'}")
            sys.exit(0)
    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --force-update / --db-version)")
    
    executor = EnteroAMRfinderPlus(cpus=args.cpus)
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        total_hits = sum(r['hit_count'] for r in results.values())
        high = sum(1 for r in results.values() for g in [h.get('gene_symbol') for h in r['hits'] if h.get('gene_symbol')] if g in executor.high_risk_genes)
        crit = sum(1 for r in results.values() for g in [h.get('gene_symbol') for h in r['hits'] if h.get('gene_symbol')] if g in executor.critical_risk_genes)
        executor.logger.info("\n📊 FINAL SUMMARY:")
        executor.logger.info(f"   Genomes: {len(results)}")
        executor.logger.info(f"   Total AMR hits: {total_hits}")
        executor.logger.info(f"   High-risk genes: {high}")
        executor.logger.info(f"   CRITICAL risk genes: {crit}")
        executor.logger.info(f"   Database: {executor.metadata['database_version']}")
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
